"""End-to-end stages. Each stage reads the data contract and writes tables under results/<run>/."""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import __version__, contracts
from .effects import two_by_two_example
from .isr_translation import dte
from .modeling import models
from .tet2_transcriptome import analysis as tet2
from .validation import evaluate


def load_cfg(path) -> dict:
    cfg = yaml.safe_load(Path(path).read_text())
    cfg["_path"] = str(path)
    cfg["out"] = Path("results") / cfg["run_name"]
    cfg["out"].mkdir(parents=True, exist_ok=True)
    return cfg


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def stage_smoke_data(cfg):
    from .data_ingest.synthetic import generate
    generate(cfg["data_root"], **cfg["synthetic"])


def stage_tet2(cfg):
    out = cfg["out"]
    r = tet2.run_all(contracts.load_tet2_studies(cfg["data_root"]), cfg["target_compartment"])
    r["per_study"].assign(evidence="observed").to_parquet(out / "tet2_gene_contrasts_per_study.parquet", index=False)
    r["meta"].assign(evidence="observed", compartment=str(cfg["target_compartment"])).to_parquet(
        out / "tet2_gene_contrasts_meta.parquet", index=False)
    r["qc"].to_csv(out / "tet2_sample_qc.tsv", sep="\t", index=False)
    r["concordance"].to_csv(out / "tet2_cross_study_spearman.csv")
    r["composition"].to_csv(out / "tet2_composition.tsv", sep="\t", index=False)
    return r


def _isr_inputs(cfg):
    root = Path(cfg["data_root"])
    ctx, smp, cnt = contracts.load_isr(root)
    genes = pd.read_parquet(root / "genes.parquet")
    programs = pd.read_csv(root / "programs.tsv", sep="\t")
    return ctx, smp, cnt, genes, programs


def stage_isr(cfg):
    out = cfg["out"]
    ctx, smp, cnt, genes, programs = _isr_inputs(cfg)
    sig = dte.build_signatures(ctx[contracts.CONTEXT_COLS], smp, cnt)
    sig.to_parquet(out / "isr_translation_signatures.parquet", index=False)
    dte.library_qc(cnt, smp).to_csv(out / "isr_library_qc.tsv", sep="\t", index=False)
    dte.pooled_signature(sig).to_parquet(out / "isr_pooled_dte.parquet", index=False)
    _specificity(cfg, sig, genes, programs)
    return sig


def _specificity(cfg, sig, genes, programs):
    """Reference-gene check for every context, and comparison of each stress signature with its
    ISR-blocked counterpart (ISRIB, eIF2a-S51A, Hri-KO). Descriptive; never used for fitting."""
    out, root = cfg["out"], Path(cfg["data_root"])
    if not (root / "isr/contexts_all.tsv").exists():
        return
    ctx = pd.read_csv(root / "isr/contexts_all.tsv", sep="\t")
    smp = pd.read_csv(root / "isr/samples_all.tsv", sep="\t")
    cnt = pd.read_parquet(root / "isr/counts_all.parquet")
    extra = ctx[~ctx.context_id.isin(sig.context_id)]
    allsig = pd.concat([sig.assign(role="train"), dte.build_signatures(
        extra, smp[smp.context_id.isin(extra.context_id)], cnt).merge(extra[["context_id"]], on="context_id")],
        ignore_index=True)
    allsig["role"] = allsig.context_id.map(ctx.set_index("context_id").role)
    allsig.to_parquet(out / "isr_translation_signatures_all.parquet", index=False)
    sym = genes.set_index("gene_id").symbol if "symbol" in genes else None
    ref = programs.loc[programs.program == cfg["reference_program"], "gene_id"]
    r = allsig[allsig.gene_id.isin(ref)].copy()
    if sym is not None:
        r.insert(0, "symbol", r.gene_id.map(sym))
    med = allsig.groupby("context_id").dTE_log2fc.median()
    r["dTE_minus_context_median"] = r.dTE_log2fc - r.context_id.map(med)
    r[["symbol", "context_id", "role", "dRNA_log2fc", "dRibo_log2fc", "dTE_log2fc", "dTE_minus_context_median",
       "dTE_se", "dTE_fdr"]].to_csv(out / "reference_set_check.csv", index=False)
    rows = []
    pairs = ctx[ctx.role == "specificity_control"]
    for _, b in pairs.iterrows():
        cand = ctx[(ctx.study_id == b.study_id) & (ctx.role != "specificity_control")]
        if cand.empty:
            continue
        a = cand.iloc[0]
        x = allsig[allsig.context_id == a.context_id].set_index("gene_id")
        y = allsig[allsig.context_id == b.context_id].set_index("gene_id")
        g = x.index.intersection(y.index)
        xr, yr = x.loc[g.intersection(ref)], y.loc[g.intersection(ref)]
        sx = x.loc[g][x.loc[g].dTE_fdr < 0.05]
        rows.append(dict(stress_context=a.context_id, stress_role=a.role, blocked_context=b.context_id,
                         n_genes=len(g), pearson_dTE_stress_vs_blocked=float(np.corrcoef(x.loc[g].dTE_log2fc, y.loc[g].dTE_log2fc)[0, 1]),
                         n_sig_dTE_stress=int((x.loc[g].dTE_fdr < 0.05).sum()), n_sig_dTE_blocked=int((y.loc[g].dTE_fdr < 0.05).sum()),
                         frac_sig_stress_genes_with_smaller_abs_dTE_when_blocked=(
                             float((np.abs(y.loc[sx.index].dTE_log2fc) < np.abs(sx.dTE_log2fc)).mean()) if len(sx) else np.nan),
                         ref_set_mean_dTE_stress=float(xr.dTE_log2fc.mean() - x.dTE_log2fc.median()),
                         ref_set_mean_dTE_blocked=float(yr.dTE_log2fc.mean() - y.dTE_log2fc.median()),
                         n_ref_genes=len(xr)))
    pd.DataFrame(rows).to_csv(out / "isr_specificity.csv", index=False)


def _context_table(cfg, sig, source=None):
    ctx, smp, cnt, genes, programs = _isr_inputs(cfg)
    st = cfg["state_programs"]
    c = ctx.set_index("context_id").join(models.context_states(cnt, smp, programs, st))
    c["intensity_metadata"] = c["intensity"]
    ref = programs.loc[programs.program == cfg["reference_program"], "gene_id"]
    eng = models.isr_engagement(sig, ref, cfg.get("engagement_min_dte", -np.inf), cfg.get("engagement_alpha", 1.01))
    c = c.join(eng)
    c.reset_index().to_csv(cfg["out"] / "isr_engagement.tsv", sep="\t", index=False)
    if "engagement_min_dte" in cfg:
        c = c[c.engaged]
        sig = sig[sig.context_id.isin(c.index)]
    c["intensity_estimated"] = models.estimate_intensity(sig, ref)
    c["intensity"] = c[f"intensity_{source or cfg['intensity_source']}"]
    c["intensity"] = c["intensity"] / c["intensity"].mean()
    return c, genes, programs, ref


def stage_baselines(cfg):
    """M0 and M1 plus the worked 2x2 identifiability example."""
    out = cfg["out"]
    ex = two_by_two_example()
    ex["cells"].to_csv(out / "example_2x2_cells.csv", index=False)
    ex["summary"].to_csv(out / "example_2x2_summary.csv", index=False)
    r = stage_tet2(cfg)
    sig = pd.read_parquet(out / "isr_translation_signatures.parquet")
    c, genes, programs, _ = _context_table(cfg, sig)
    models.m0_programs(r["expr"], programs).to_csv(out / "m0_program_contrasts.csv", index=False)
    s2 = sig.drop(columns=["intensity"]).merge(c[["intensity"]], left_on="context_id", right_index=True)
    pooled = dte.pooled_signature(s2, per_unit_intensity=True)
    g, p = models.m1_predict(r["expr"], pooled, programs, cfg["isr_intensities"], cfg["n_mc"], cfg["seed"])
    g.to_parquet(out / "m1_gene_predictions.parquet", index=False)
    p.to_csv(out / "m1_program_predictions.csv", index=False)
    # sensitivity of M1 program predictions to source choice: leave one ISR study out
    rows = []
    for s in c.study_id.unique():
        keep = c.index[c.study_id != s]
        pl = dte.pooled_signature(s2[s2.context_id.isin(keep)], per_unit_intensity=True)
        _, pp = models.m1_predict(r["expr"], pl, programs, {"moderate": cfg["isr_intensities"]["moderate"]},
                                  n_mc=50, seed=cfg["seed"])
        rows.append(pp.assign(isr_study_left_out=s))
    pd.concat(rows).to_csv(out / "m1_source_sensitivity.csv", index=False)


def stage_evaluate(cfg):
    """M2/M3, held-out evaluation, ablations, leakage, genotype prediction, direct-holdout check."""
    out = cfg["out"]
    sig = pd.read_parquet(out / "isr_translation_signatures.parquet")
    r = tet2.run_all(contracts.load_tet2_studies(cfg["data_root"]), cfg["target_compartment"])
    st = cfg["state_programs"]
    sig_all = sig
    c, genes, programs, ref = _context_table(cfg, sig)
    sig = sig[sig.context_id.isin(c.index)]
    c.reset_index().to_csv(out / "isr_context_table.tsv", sep="\t", index=False)
    evaluate.unengaged_check(sig_all, sig, genes, c, st, programs).to_csv(out / "unengaged_context_check.csv", index=False)
    ev = evaluate.run_heldout(sig, genes, c, st, programs, ref, cfg["seed"])
    ev.to_csv(out / "heldout_metrics.csv", index=False)
    evaluate.gene_family_leakage(sig, genes, c, seed=cfg["seed"]).to_csv(out / "gene_family_leakage.csv", index=False)
    evaluate.factor_stability(sig, genes, c, st, seed=cfg["seed"]).to_csv(out / "m3_factor_stability.csv", index=False)
    preds, hold_rows = [], []
    hold = Path(cfg["data_root"]) / "tet2_isr_holdout"
    obs = None
    if hold.exists():
        hs = pd.read_csv(hold / "samples.tsv", sep="\t")
        hc = contracts.check_counts(pd.read_parquet(hold / "counts.parquet"), hs, "holdout")
        obs = dte.genotype_by_isr_dte(hc, hs)
        obs.to_parquet(out / "holdout_observed_genotype_x_isr.parquet", index=False)
    m1 = pd.read_parquet(out / "m1_gene_predictions.parquet")
    need = cfg.get("min_contexts_for_state_model", 0)
    if len(c) < need:
        msg = (f"M2/M3 genotype prediction NOT RUN: {len(c)} ISR-engaged contexts from {c.study_id.nunique()} studies; "
               f"the state-dependent model has {len(st)} state covariates and needs >= {need} contexts "
               "(config: min_contexts_for_state_model). Model-based I_TE is not identifiable from these data.\n")
        (out / "m2_genotype_predictions.NOT_SUPPORTED.txt").write_text(msg)
        (out / "holdout_validation.NOT_AVAILABLE.txt").write_text(
            "No direct genotype x ISR translatome exists in public data (GEO search 2026-11-12: 0 series).\n")
        print("[tet2isr] " + msg.strip())
        return
    for src in ("estimated", "metadata"):
        cs, _, _, _ = _context_table(cfg, sig_all, src)
        if cs.intensity.isna().any():
            continue
        p, states, _ = models.m2_genotype_prediction(sig, genes, cs, st, r["expr"], programs,
                                                     cfg["isr_intensities"], cfg["n_boot"], cfg["seed"])
        preds.append(p.assign(intensity_source=src))
        states.reset_index(names="genotype").to_csv(out / "tet2_state_scores.tsv", sep="\t", index=False)
        if obs is not None:
            hold_rows.append(evaluate.genotype_holdout(obs, m1, p).assign(intensity_source=src))
    allp = pd.concat(preds)
    allp.to_parquet(out / "m2_genotype_predictions.parquet", index=False)
    if obs is not None:
        evaluate.holdout_program_level(obs, allp, programs).to_csv(out / "holdout_program_level.csv", index=False)
    if hold_rows:
        pd.concat(hold_rows).to_csv(out / "holdout_validation.csv", index=False)
    else:
        (out / "holdout_validation.NOT_AVAILABLE.txt").write_text(
            "No direct genotype x ISR translatome in this data root. I_TE predictions are unvalidated.\n")


def stage_intersect(cfg):
    """TET2 transcriptome programs x stress translatomes (needs analyze-isr, fit-baselines, evaluate)."""
    from . import intersection
    out, root = cfg["out"], Path(cfg["data_root"])
    if not (out / "isr_translation_signatures_all.parquet").exists():
        return
    species = pd.read_csv(root / "isr/contexts_all.tsv", sep="\t").set_index("study_id").species.to_dict()
    for d in sorted((root / "tet2").iterdir()):
        sp = pd.read_csv(d / "samples.tsv", sep="\t").species.iloc[0]
        species[d.name] = "mouse" if "mouse" in sp else "human"
    r = intersection.run(out, pd.read_parquet(root / "genes.parquet"), species, seed=cfg["seed"])
    (out / "figures").mkdir(exist_ok=True)
    intersection.figure(r["pairs"], out / "figures/fig9_tet2_x_translatome.png")
    (out / "intersection_null_q95.txt").write_text(f"{r['null_q95']:.4f}\n")


def write_provenance(cfg):
    out = cfg["out"]
    files = {str(p): _sha(p) for p in sorted(out.rglob("*")) if p.is_file() and p.name != "provenance.json"}
    inputs = {str(p): _sha(p) for p in sorted(Path(cfg["data_root"]).rglob("*")) if p.is_file()}
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except OSError:
        commit = ""
    import matplotlib, scipy  # noqa
    prov = dict(run_name=cfg["run_name"], generated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                config={k: (str(v) if isinstance(v, Path) else v) for k, v in cfg.items()},
                data_is_synthetic=bool(cfg.get("synthetic")), git_commit=commit,
                software=dict(tet2isr=__version__, python=sys.version.split()[0], platform=platform.platform(),
                              numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                              matplotlib=matplotlib.__version__),
                input_sha256=inputs, output_sha256=files)
    (out / "provenance.json").write_text(json.dumps(prov, indent=1))
