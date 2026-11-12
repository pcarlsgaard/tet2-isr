"""Figures, priority table and HTML report. Everything is regenerated from results tables."""
from __future__ import annotations

import base64
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

WT, KO, ACC, GREY, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#b9b8b2", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlesize": 10, "axes.titleweight": "bold",
                     "axes.titlelocation": "left", "figure.facecolor": "white", "legend.frameon": False,
                     "savefig.dpi": 200, "savefig.bbox": "tight", "lines.linewidth": 2})


def _save(fig, path):
    fig.savefig(path)
    plt.close(fig)
    return path


def fig1_evidence_map(manifest: pd.DataFrame, path):
    caps = ["supports_baseline_tet2", "supports_isr_transcription", "supports_isr_translation",
            "supports_direct_genotype_x_isr"]
    lab = ["basal TET2\ncontrast", "ISR\ntranscription", "ISR translation\n(paired RNA+Ribo)", "direct\nTET2 x ISR TE"]
    m = manifest[manifest.accession.notna() & (manifest.source_class != "literature_prior")].reset_index(drop=True)
    col = {"yes": WT, "no": "#ecebe7", "unverified": "#f3d9a4", "claimed_unverified": "#f3d9a4"}
    fig, ax = plt.subplots(figsize=(7.2, 0.36 * len(m) + 1.3))
    for i, r in m.iterrows():
        for j, c in enumerate(caps):
            v = str(r[c])
            ax.add_patch(plt.Rectangle((j + 0.04, i + 0.08), 0.92, 0.84, color=col.get(v, "#ecebe7")))
            ax.text(j + 0.5, i + 0.5, {"claimed_unverified": "claimed,\nunverified"}.get(v, v), ha="center",
                    va="center", fontsize=7, color="white" if v == "yes" else INK)
    ax.set_xlim(0, 4); ax.set_ylim(len(m), 0)
    ax.set_xticks(np.arange(4) + 0.5, lab); ax.xaxis.tick_top(); ax.tick_params(length=0)
    ax.set_yticks(np.arange(len(m)) + 0.5, [f"{a}  ({s})" for a, s in zip(m.accession, m.species)])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Evidence map: what each candidate public dataset can support\n"
                 "(amber = claimed in the brief, sample-level metadata not verifiable this run)", pad=34)
    return _save(fig, path)


def fig2_basal(m0: pd.DataFrame, path):
    progs = list(m0.groupby("program").log2fc.mean().sort_values().index)
    fig, ax = plt.subplots(figsize=(6, 0.34 * len(progs) + 1.2))
    studies = list(m0.study_id.unique())
    for k, (s, c) in enumerate(zip(studies, (WT, ACC, KO))):
        d = m0[m0.study_id == s].set_index("program").loc[progs]
        y = np.arange(len(progs)) + (k - (len(studies) - 1) / 2) * 0.22
        ax.errorbar(d.log2fc, y, xerr=[d.log2fc - d.ci_lo, d.ci_hi - d.log2fc], fmt="o", color=c, ms=5,
                    lw=1.5, label=f"{s} (n={d.n_WT.iloc[0]} vs {d.n_TET2.iloc[0]})")
    ax.axvline(0, color=GREY, lw=1)
    ax.set_yticks(range(len(progs)), progs)
    ax.set_xlabel("program score, log2(TET2 / WT), 95% CI over biological replicates")
    ax.set_title("M0: basal RNA programs by genotype (HSC compartment)")
    ax.legend(loc="lower right", fontsize=7)
    return _save(fig, path)


def fig3_rna_vs_te(sig, qc, ref, path):
    cid = sig.groupby("context_id").intensity.first().idxmax()
    d = sig[sig.context_id == cid]
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.8), gridspec_kw=dict(width_ratios=[1, 1.25]))
    r = d.gene_id.isin(ref)
    a.scatter(d.dRNA_log2fc[~r], d.dTE_log2fc[~r], s=5, color=GREY, label="other genes", rasterized=True)
    a.scatter(d.dRNA_log2fc[r], d.dTE_log2fc[r], s=12, color=KO, label="reference uORF targets")
    a.axhline(0, color=MUTED, lw=0.6); a.axvline(0, color=MUTED, lw=0.6)
    a.set_xlabel("RNA change, log2(stress / control)"); a.set_ylabel("dTE, log2")
    a.set_title(f"RNA vs TE change\n{cid}"); a.legend(fontsize=7)
    q = qc.pivot(index="context_id", columns="assay", values="min_within_condition_r")
    y = np.arange(len(q))
    b.scatter(q.rna, y - 0.12, color=WT, s=28, label="RNA-seq"); b.scatter(q.ribo, y + 0.12, color=KO, s=28, label="Ribo-seq")
    b.set_yticks(y, [c.split("|", 1)[1] for c in q.index], fontsize=7)
    b.set_xlabel("minimum replicate correlation within condition (log2 CPM)")
    b.set_title("Count-level QC (read-level QC not available\nfrom processed counts)"); b.legend(fontsize=7)
    fig.tight_layout()
    return _save(fig, path)


def fig4_obs_pred(sig, ctx, genes, state_cols, path):
    from ..modeling.models import make_model
    study = ctx.groupby("study_id").size().idxmax()
    te = ctx[ctx.study_id == study]
    tr = ctx[ctx.study_id != study]
    m = make_model("M2_+cell_state", genes, state_cols).fit(sig[sig.context_id.isin(tr.index)], genes, tr)
    P = m.predict(genes, te)
    fig, axes = plt.subplots(1, len(te), figsize=(4.2 * len(te), 3.9), squeeze=False)
    for ax, c in zip(axes[0], te.index):
        o = sig[sig.context_id == c].set_index("gene_id").dTE_log2fc
        p = P[c].reindex(o.index)
        o, p = o - o.median(), p - p.median()
        ax.scatter(p, o, s=5, color=WT, alpha=0.5, rasterized=True)
        lim = [min(p.min(), o.min()), max(p.max(), o.max())]
        ax.plot(lim, lim, color=MUTED, lw=1, ls="--")
        ax.set_xlabel("predicted dTE (study held out, dose-blind)"); ax.set_ylabel("observed dTE")
        ax.set_title(f"{c.split('|', 1)[1]}\nPearson r = {np.corrcoef(o, p)[0, 1]:.2f}")
    fig.suptitle(f"Held-out study: {study} (never seen in training)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout()
    return _save(fig, path)


def fig5_ablation(ev, path):
    order = ["NULL_permuted_labels", "B0_transcriptome_only/uniform_shift", "B2_reference_gene_set",
             "B1_nearest_context", "M2a_sequence_only", "M1_gene_signature", "M1b_signature+dose",
             "M2_+cell_state", "M2c_+stress_class", "M2d_+epigenome", "M3_low_rank"]
    splits = ["held_out_stressor", "held_out_cell_type", "held_out_study"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), sharey=True)
    for ax, met, lab in zip(axes, ("pearson", "rmse"), ("Pearson r (held-out contexts)", "RMSE, log2 (dose-blind)")):
        for k, (sp, c) in enumerate(zip(splits, (WT, KO, ACC))):
            d = ev[ev.split == sp].groupby("model")[met].agg(["mean", "min", "max"]).reindex(order)
            y = np.arange(len(order)) + (k - 1) * 0.22
            ax.errorbar(d["mean"], y, xerr=[d["mean"] - d["min"], d["max"] - d["mean"]], fmt="o", ms=4.5,
                        lw=1.2, color=c, label=sp.replace("_", " "))
        ax.set_xlabel(lab + "\nmean and range over held-out contexts")
        ax.axvline(0, color=GREY, lw=0.8)
    axes[0].set_yticks(range(len(order)), order)
    axes[0].set_title("Ablation ladder and baselines")
    axes[0].set_title("Ablation ladder and baselines (Pearson r is undefined for the constant B0 prediction)")
    h, l = axes[1].get_legend_handles_labels()
    fig.tight_layout()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=8, bbox_to_anchor=(0.5, -0.06))
    return _save(fig, path)


def fig6_genotype(m1p, hp, path):
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4), sharey=True)
    d = m1p[m1p.isr_intensity == "moderate"]
    progs = sorted(d.program.unique())
    for k, (s, c) in enumerate(zip(d.tet2_study.unique(), (WT, ACC))):
        x = d[d.tet2_study == s].set_index("program").loc[progs]
        y = np.arange(len(progs)) + (k - 0.5) * 0.25
        a.errorbar(x.I_S_log2_program, y, xerr=[x.I_S_log2_program - x.I_S_log2_program_lo,
                                                x.I_S_log2_program_hi - x.I_S_log2_program],
                   fmt="o", ms=4.5, lw=1.2, color=c, label=s)
    a.axvline(0, color=GREY, lw=1); a.set_yticks(range(len(progs)), progs)
    a.set_xlabel("program-level I_S, log2 (95% MC interval)")
    a.set_title("M1: substrate-only interaction\n(gene-level I_TE = 0 by assumption)"); a.legend(fontsize=7, loc="lower right")
    if hp is not None and len(hp):
        h = hp[hp.tet2_study == hp.tet2_study.iloc[0]]
        for k, (src, c) in enumerate((("estimated", KO), ("metadata", ACC))):
            x = h[h.intensity_source == src].set_index("program").reindex(progs)
            y = np.arange(len(progs)) + (k - 0.5) * 0.25
            b.errorbar(x.predicted_I_TE_M2, y, xerr=[x.predicted_I_TE_M2 - x.predicted_lo,
                                                     x.predicted_hi - x.predicted_I_TE_M2],
                       fmt="o", ms=4.5, lw=1.2, color=c, label=f"M2 prediction, dose = {src}")
        x = h[h.intensity_source == "estimated"].set_index("program").reindex(progs)
        b.errorbar(x.observed_I_TE, np.arange(len(progs)), xerr=1.96 * x.observed_se, fmt="D", ms=4,
                   lw=1.2, color=INK, label="direct measurement (held-out 2x2 translatome)")
    b.axvline(0, color=GREY, lw=1)
    b.set_xlabel("program-mean I_TE, log2 (95% interval)")
    b.set_title("M2: model-based I_TE vs direct measurement")
    fig.tight_layout()
    h, l = b.get_legend_handles_labels()
    if h:
        fig.legend(h, l, loc="lower right", ncol=3, fontsize=8, bbox_to_anchor=(0.99, -0.07))
    return _save(fig, path)


def fig7_hypotheses(path):
    x = np.linspace(0, 3, 200)
    hill = lambda x, ec, top: top * x ** 2 / (ec ** 2 + x ** 2)
    spec = [("H0  no genotype-specific TE difference", 1.0, 1.0), ("H1  enhanced adaptive ISR", 1.0, 1.45),
            ("H2  blunted adaptive ISR", 1.0, 0.6), ("H3  shifted dose threshold", 0.55, 1.0)]
    fig, axes = plt.subplots(1, 4, figsize=(11, 2.9), sharey=True)
    for ax, (t, ec, top) in zip(axes, spec):
        ax.plot(x, hill(x, 1.0, 1.0), color=WT, label="WT")
        ax.plot(x, hill(x, ec, top), color=KO, ls="--", label="TET2-deficient")
        ax.set_title(t, fontsize=8.5); ax.set_xlabel("ISR input (p-eIF2α, matched across genotypes)")
        ax.set_xticks([]); ax.set_yticks([])
    axes[0].set_ylabel("dTE of uORF-target set"); axes[0].legend(fontsize=7)
    fig.suptitle("Conceptual sketch, not data: the four competing hypotheses differ in dose-response SHAPE, "
                 "so a dose series at matched p-eIF2α separates them", x=0.01, ha="left", fontsize=9)
    fig.tight_layout()
    return _save(fig, path)


# ---------------------------------------------------------------- priority tables
def priority_table_run(out: Path, synthetic: bool) -> pd.DataFrame:
    """Model-derived priority table for the current run (program level + top genes)."""
    m1p = pd.read_csv(out / "m1_program_predictions.csv")
    m2 = pd.read_parquet(out / "m2_genotype_predictions.parquet")
    programs = pd.read_csv(Path(out / "_programs.tsv"), sep="\t")
    src = "SYNTHETIC smoke-test data" if synthetic else "see provenance.json"
    rows = []
    m2 = m2.merge(programs, on="gene_id")
    agg = (m2.groupby(["program", "isr_intensity", "intensity_source"])
           .agg(dTE_WT=("dTE_WT", "mean"), dTE_TET2=("dTE_TET2", "mean"), I_TE=("I_TE", "mean"),
                lo=("I_TE_lo", "mean"), hi=("I_TE_hi", "mean"), sign=("sign_stability", "mean"),
                rank=("rank_stability_spearman", "mean"), extrap=("extrapolation_flag", "max")).reset_index())
    s1 = m1p.groupby(["program", "isr_intensity"])[["I_S_log2_program", "I_S_log2_program_lo",
                                                     "I_S_log2_program_hi"]].mean().reset_index()
    for _, r in agg.merge(s1, on=["program", "isr_intensity"]).iterrows():
        rows.append({
            "gene_or_program": r.program, "isr_dose": r.isr_intensity, "time_phase": "early (translational)",
            "intensity_source": r.intensity_source,
            "predicted_direction_WT": "up" if r.dTE_WT > 0.1 else "down" if r.dTE_WT < -0.1 else "flat",
            "predicted_direction_TET2": "up" if r.dTE_TET2 > 0.1 else "down" if r.dTE_TET2 < -0.1 else "flat",
            "predicted_I_TE_M2": round(r.I_TE, 3), "I_TE_interval": f"[{r.lo:.3f}, {r.hi:.3f}]",
            "I_TE_defensible": bool((r.lo > 0 or r.hi < 0) and not r.extrap),
            "predicted_I_S_log2_M1_substrate_only": round(r.I_S_log2_program, 3),
            "I_S_interval": f"[{r.I_S_log2_program_lo:.3f}, {r.I_S_log2_program_hi:.3f}]",
            "evidence_category": "transferred", "source_ids": src,
            "model_stability": f"sign {r.sign:.2f}; rank rho {r['rank']:.2f}; extrapolation={bool(r.extrap)}",
            "plausible_mechanism": "state-dependent uORF/TOP response learned across contexts (M2)",
            "alternative_explanation": "dose mis-calibration; substrate-only effect (M1); collinear state programs",
            "recommended_discriminating_readout": "set-level I_TE from 2x2 paired RNA/Ribo at matched p-eIF2a"})
    t = pd.DataFrame(rows).sort_values("predicted_I_TE_M2", key=np.abs, ascending=False)
    t.to_csv(out / "prediction_priority_table.csv", index=False)
    return t


def priority_table_public(prior_path="configs/prior_hypotheses.tsv", out="results/prediction_priority_table.csv"):
    t = pd.read_csv(prior_path, sep="\t")
    t.insert(0, "status", "PRIOR-ONLY: no public-data model fit (data hosts blocked in the build session); "
                          "no numeric effect is estimated")
    Path(out).parent.mkdir(exist_ok=True)
    t.to_csv(out, index=False)
    return t


# ---------------------------------------------------------------- HTML
def _img(p):
    return f'<img src="data:image/png;base64,{base64.b64encode(Path(p).read_bytes()).decode()}" style="max-width:100%">'


def _tbl(df, n=None, r=3):
    d = df.head(n) if n else df
    return d.round(r).to_html(index=False, border=0, classes="t", na_rep="–")


CSS = """body{font:15px/1.5 system-ui,sans-serif;color:#0b0b0b;max-width:1060px;margin:2em auto;padding:0 1em}
h1{font-size:1.6em}h2{margin-top:2em;border-bottom:1px solid #ddd;padding-bottom:.2em}
.banner{background:#fff3cd;border:1px solid #e0b94a;padding:.8em 1em;border-radius:6px;font-weight:600}
.t{border-collapse:collapse;font-size:12px;margin:.6em 0;display:block;overflow-x:auto}
.t th,.t td{padding:3px 8px;border-bottom:1px solid #e6e6e6;text-align:left;white-space:nowrap}
.t th{background:#f5f5f3}.note{color:#52514e;font-size:13px}code{background:#f3f3f1;padding:1px 4px;border-radius:3px}"""


def build(cfg) -> Path:
    out = Path(cfg["out"])
    figs = out / "figures"
    figs.mkdir(exist_ok=True)
    synthetic = bool(cfg.get("synthetic"))
    root = Path(cfg["data_root"])
    programs = pd.read_csv(root / "programs.tsv", sep="\t")
    programs.to_csv(out / "_programs.tsv", sep="\t", index=False)
    genes = pd.read_parquet(root / "genes.parquet")
    manifest = pd.read_csv("data/dataset_manifest.tsv", sep="\t")
    sig = pd.read_parquet(out / "isr_translation_signatures.parquet")
    ctx = pd.read_csv(out / "isr_context_table.tsv", sep="\t").set_index("context_id")
    qc = pd.read_csv(out / "isr_library_qc.tsv", sep="\t")
    m0 = pd.read_csv(out / "m0_program_contrasts.csv")
    m1p = pd.read_csv(out / "m1_program_predictions.csv")
    ev = pd.read_csv(out / "heldout_metrics.csv")
    hp = pd.read_csv(out / "holdout_program_level.csv") if (out / "holdout_program_level.csv").exists() else None
    hv = pd.read_csv(out / "holdout_validation.csv") if (out / "holdout_validation.csv").exists() else None
    ref = set(programs.loc[programs.program == cfg["reference_program"], "gene_id"])
    f = [fig1_evidence_map(manifest, figs / "fig1_evidence_map.png"),
         fig2_basal(m0, figs / "fig2_basal_programs.png"),
         fig3_rna_vs_te(sig, qc, ref, figs / "fig3_rna_vs_te_qc.png"),
         fig4_obs_pred(sig, ctx, genes, cfg["state_programs"], figs / "fig4_heldout_obs_vs_pred.png"),
         fig5_ablation(ev, figs / "fig5_ablation.png"),
         fig6_genotype(m1p, hp, figs / "fig6_genotype_by_isr.png"),
         fig7_hypotheses(figs / "fig7_hypotheses.png")]
    pri = priority_table_run(out, synthetic)
    pub = priority_table_public()
    summ = (ev.groupby(["split", "model"])[["pearson", "spearman", "rmse", "sign_concordance_sig",
                                            "top50_precision", "program_recovery_r", "coverage_50",
                                            "coverage_80", "coverage_95"]].mean().reset_index())
    summ.to_csv(out / "heldout_summary.csv", index=False)
    seen = (ev[(ev.split == "held_out_stressor") & ev.model.isin(["M1b_signature+dose", "M2_+cell_state"])]
            .groupby(["model", "stressor_class", "cell_type_in_training"])[["pearson", "rmse"]].mean().reset_index())
    banner = ("SYNTHETIC DEMONSTRATION. Every number, figure and gene ID (SYN…) below comes from simulated data "
              "with planted effects. It shows that the software works and what the design can and cannot "
              "identify. It says nothing about real TET2 biology." if synthetic else
              "Public-data run. Predictions are transfer-based unless marked observed.")
    h = [f"<html><head><meta charset='utf-8'><title>TET2 x ISR report ({cfg['run_name']})</title><style>{CSS}</style></head><body>",
         f"<h1>TET2 × ISR translatomic transfer — research report ({cfg['run_name']})</h1>",
         f"<p class='banner'>{banner}</p>",
         "<p class='note'>Regenerate with <code>make smoke-test</code>. Methods: docs/methods.md. Limits: "
         "docs/assumptions_and_limits.md. Inputs and checksums: provenance.json in the run directory.</p>",
         "<h2>1. Datasets and evidence map</h2>", _img(f[0]),
         "<p class='note'>Status column of data/dataset_manifest.tsv records how each field was verified. "
         "No public source is known to provide a direct TET2 × ISR translatome.</p>",
         _tbl(manifest[["accession", "species", "assay", "verification_status", "relevance_rank", "confounding_risk"]]),
         "<h2>2. Identifiability: the 2×2 example</h2>",
         "<p>Gene A doubles in mRNA with TET2 loss; ISR doubles its TE in both genotypes. Linear output differs "
         "(I_S_linear = 10 units) while I_TE = 0 and gene-level log I_S = 0. The two-gene pathway shows a non-zero "
         "<em>aggregate</em> log interaction purely because the abundance weights differ.</p>",
         _tbl(pd.read_csv(out / "example_2x2_summary.csv"), r=4),
         "<h2>3. Basal genotype contrasts (M0)</h2>", _img(f[1]),
         _tbl(pd.read_csv(out / "tet2_sample_qc.tsv", sep="\t")),
         "<h2>4. ISR translation signatures and QC</h2>", _img(f[2]),
         "<h2>5. Held-out prediction</h2>", _img(f[3]), _img(f[4]),
         "<p class='note'>Mean over held-out contexts. coverage_XX = fraction of held-out genes inside the "
         "nominal XX% predictive interval.</p>", _tbl(summ),
         "<p class='note'>Held-out stressor, split by stressor class and by whether the test cell type was seen:</p>",
         _tbl(seen),
         "<p class='note'>Gene-family leakage check (sequence-only model, unseen genes):</p>",
         _tbl(pd.read_csv(out / "gene_family_leakage.csv")),
         "<p class='note'>M3 factor stability (delete-2-contexts jackknife):</p>",
         _tbl(pd.read_csv(out / "m3_factor_stability.csv")),
         "<h2>6. Genotype × ISR predictions</h2>", _img(f[5])]
    if hv is not None:
        h += ["<p class='note'>Gene-level comparison with the held-out direct 2×2 translatome (moderate dose):</p>",
              _tbl(hv[hv.isr_intensity == "moderate"][["intensity_source", "model", "tet2_study", "target",
                                                        "pearson", "rmse", "n_sig_obs", "top50_precision",
                                                        "coverage_95"]]),
              "<p class='note'>Program level:</p>", _tbl(hp)]
    else:
        h.append("<p class='banner'>No direct genotype × ISR translatome is available: I_TE is an unvalidated model-based prediction.</p>")
    h += ["<h3>Run-specific priority table (top rows)</h3>", _tbl(pri, 16),
          "<h2>7. Competing hypotheses and discriminating assays</h2>", _img(f[6]),
          "<p class='note'>Real-gene hypothesis table (prior knowledge only, no fitted magnitudes): "
          "results/prediction_priority_table.csv</p>",
          _tbl(pub[["gene_or_program", "predicted_direction_WT", "predicted_direction_TET2", "predicted_I_TE",
                    "evidence_category", "recommended_discriminating_readout"]]),
          "</body></html>"]
    Path("reports").mkdir(exist_ok=True)
    name = "research_report.html" if cfg["run_name"] in ("smoke", "public") else f"research_report_{cfg['run_name']}.html"
    p = Path("reports") / name
    p.write_text("\n".join(h))
    return p
