"""Intersection of TET2-mutant transcriptome programs with stress translatomes.

Question: are the genes that TET2 loss changes at the RNA level the same genes whose translation the
ISR changes?  Parsed several ways: (i) every TET2 contrast x every translatome context, flagged
species-matched or not; (ii) TET2 'core' programs (cross-species, mouse-only, within-study
multi-population, human) against the engaged-ISR core; (iii) gene-level overlap lists.

Guards against false positives:
  * partial rank correlation removes shared technical axes (expression level, CDS length, GC, 3'UTR length);
  * ISR-blocked contexts (ISRIB, eIF2a-S51A, Hri-KO) and gate-failed stress contexts give an empirical null
    for how large a TET2 x translatome correlation gets without an ISR;
  * set tests use expression-matched random gene sets.
Everything here is association between separate datasets. None of it measures translation in TET2 cells.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

from .stats import bh_fdr, random_effects_meta

COVARS = ["seq_log_cds_len", "seq_cds_gc", "seq_log_utr3_len", "seq_utr5_gc"]


def _partial_spearman(x: pd.Series, y: pd.Series, cov: pd.DataFrame) -> tuple[float, float, int]:
    d = pd.concat([x.rename("x"), y.rename("y"), cov], axis=1, join="inner").dropna()
    if len(d) < 500:
        return np.nan, np.nan, len(d)
    r = d.rank()
    Z = np.column_stack([np.ones(len(r)), r.drop(columns=["x", "y"]).values])
    res = lambda v: v - Z @ np.linalg.lstsq(Z, v, rcond=None)[0]
    return float(sps.spearmanr(d.x, d.y)[0]), float(np.corrcoef(res(r.x.values), res(r.y.values))[0, 1]), len(d)


def tet2_cores(per: pd.DataFrame, mouse_studies, human_studies) -> pd.DataFrame:
    """Long table of TET2 'core' definitions: core, gene_id, log2fc, fdr, in_core, n_units."""
    per = per.assign(unit=per.study_id + "|" + per.compartment)

    def meta(units, name, within_study=False):
        d = per[per.unit.isin(units)]
        if within_study:   # populations of one study share animals: average first, no pooling of evidence
            e = d.groupby("gene_id").agg(log2fc=("log2fc", "mean"), se=("se", "mean"), k=("unit", "nunique"),
                                         pos=("log2fc", lambda v: (v > 0).sum()))
            e = e[e.k == len(units)]
            z = e.log2fc / e.se
            e["fdr"] = bh_fdr(2 * sps.norm.sf(np.abs(z)))
            e["in_core"] = (e.fdr < 0.05) & ((e.pos == e.k) | (e.pos == 0))
            return e.reset_index().assign(core=name, n_units=len(units))[["core", "gene_id", "log2fc", "fdr", "in_core", "n_units"]]
        est = d.pivot(index="gene_id", columns="unit", values="log2fc")
        se = d.pivot(index="gene_id", columns="unit", values="se")
        m = random_effects_meta(est.values, se.values)
        same = (np.sign(est.values) == np.sign(m.estimate.values)[:, None]).sum(axis=1)
        w = 1 / se.values ** 2
        idx = np.nanargmax(np.nan_to_num(w, nan=-1), axis=1)
        e2, s2 = est.values.copy(), se.values.copy()
        e2[np.arange(len(e2)), idx] = np.nan
        s2[np.arange(len(s2)), idx] = np.nan
        lo = random_effects_meta(e2, s2)
        core = ((m.fdr < 0.05) & (m.k == len(units)) & (same == len(units))
                & (np.sign(lo.estimate) == np.sign(m.estimate)) & (lo.pvalue < 0.05)).values
        return pd.DataFrame(dict(core=name, gene_id=est.index, log2fc=m.estimate.values, fdr=m.fdr.values,
                                 in_core=core, n_units=len(units)))

    tgt = per[per.is_target]
    out = [meta(list(tgt.unit.unique()), "cross_species_target_populations")]
    mt = list(tgt[tgt.study_id.isin(mouse_studies)].unit.unique())
    if len(mt) >= 3:
        out.append(meta(mt, "mouse_HSPC_target_populations"))
    for sid in list(mouse_studies) + list(human_studies):
        u = sorted(per[(per.study_id == sid) & ~per.compartment.str.contains("LPS|old")].unit.unique())
        if len(u) >= 2:
            out.append(meta(u, f"{sid}_all_populations(within-study)", within_study=True))
    return pd.concat(out, ignore_index=True)


def run(out: Path, genes: pd.DataFrame, species: dict, n_perm: int = 2000, seed: int = 0) -> dict:
    """species: study_id -> 'mouse' | 'human' for both TET2 and ISR studies."""
    rng = np.random.default_rng(seed)
    out = Path(out)
    sym = genes.set_index("gene_id").symbol
    cov = genes.set_index("gene_id")[COVARS]
    per = pd.read_parquet(out / "tet2_gene_contrasts_per_study.parquet")
    per["unit"] = per.study_id + "|" + per.compartment
    sig = pd.read_parquet(out / "isr_translation_signatures_all.parquet")
    eng = pd.read_csv(out / "isr_engagement.tsv", sep="\t").set_index("context_id").engaged
    cls = lambda r: ("ISR-blocked control" if r.role == "specificity_control" else
                     "low depth" if r.role == "excluded_low_depth" else
                     "ISR engaged" if eng.get(r.context_id, False) else "stress, ISR gate not passed")
    cmeta = sig.drop_duplicates("context_id").assign(context_class=lambda d: d.apply(cls, axis=1)) \
        .set_index("context_id")[["study_id", "context_class"]]
    D = sig.pivot(index="gene_id", columns="context_id", values="dTE_log2fc")
    D = D - D.median()

    # ---- (i) all pairs
    rows = []
    for u, d in per.groupby("unit"):
        x = d.set_index("gene_id").log2fc
        c2 = cov.join(d.set_index("gene_id").mean_log2cpm_WT)
        for c in D.columns:
            r, rp, n = _partial_spearman(x, D[c], c2)
            rows.append(dict(tet2_unit=u, tet2_species=species[d.study_id.iloc[0]], isr_context=c,
                             isr_species=species[cmeta.study_id[c]], context_class=cmeta.context_class[c],
                             n_genes=n, spearman=r, partial_spearman=rp))
    pairs = pd.DataFrame(rows)
    pairs["species_matched"] = pairs.tet2_species == pairs.isr_species
    null = pairs[pairs.context_class.isin(["ISR-blocked control", "stress, ISR gate not passed"])].partial_spearman.abs()
    q95 = float(null.quantile(0.95))
    pairs["exceeds_null_95"] = pairs.partial_spearman.abs() > q95
    pairs.to_csv(out / "intersection_pairwise.csv", index=False)

    # ---- (ii) core programs vs ISR axes
    mouse = [s for s in per.study_id.unique() if species[s] == "mouse"]
    human = [s for s in per.study_id.unique() if species[s] == "human"]
    cores = tet2_cores(per, mouse, human)
    cores.insert(2, "symbol", cores.gene_id.map(sym))
    engaged = [c for c in D.columns if cmeta.context_class[c] == "ISR engaged"]
    axes = {"engaged_ISR_core(mean of engaged contexts)": D[engaged].dropna().mean(axis=1)}
    axes.update({c: D[c].dropna() for c in D.columns})
    expr = per.groupby("gene_id").mean_log2cpm_WT.mean()
    rows = []
    for core, d in cores.groupby("core"):
        for direction, sel in (("up_in_TET2", d[d.in_core & (d.log2fc > 0)]), ("down_in_TET2", d[d.in_core & (d.log2fc < 0)])):
            for ax_name, ax in axes.items():
                bg = ax.index.intersection(d.gene_id)
                g = bg.intersection(sel.gene_id)
                row = dict(tet2_core=core, direction=direction, n_core_genes=len(sel), isr_axis=ax_name,
                           context_class=cmeta.context_class.get(ax_name, "ISR engaged"), n_measured=len(g))
                if len(g) >= 10:
                    # expression-matched permutation null for the set mean
                    dec = pd.qcut(expr.reindex(bg).rank(method="first"), 10, labels=False)
                    pools = {k: np.where(dec.values == k)[0] for k in range(10)}
                    need = dec.loc[g].value_counts()
                    v = ax.loc[bg].values
                    sims = np.array([np.concatenate([v[rng.choice(pools[k], n_k, replace=False)] for k, n_k in need.items()]).mean()
                                     for _ in range(n_perm)])
                    obs = float(ax.loc[g].mean())
                    row.update(mean_dTE=obs, null_mean=float(sims.mean()), z=float((obs - sims.mean()) / sims.std()),
                               p_perm=float((np.sum(np.abs(sims - sims.mean()) >= abs(obs - sims.mean())) + 1) / (n_perm + 1)))
                rows.append(row)
    sets = pd.DataFrame(rows)
    sets["fdr"] = bh_fdr(sets.p_perm.values) if "p_perm" in sets else np.nan
    sets.to_csv(out / "intersection_core_sets.csv", index=False)
    cores.to_csv(out / "tet2_core_programs.csv", index=False)

    # ---- (iii) gene-level overlap: TET2 core genes x ISR translation core
    Y = D[engaged]
    isr = pd.DataFrame(dict(isr_dTE_mean=Y.mean(axis=1), isr_n_contexts=Y.notna().sum(axis=1),
                            isr_min_abs=Y.abs().min(axis=1),
                            isr_consistent=(Y.notna().all(axis=1)) & ((Y > 0).all(axis=1) | (Y < 0).all(axis=1))))
    isr["isr_core"] = isr.isr_consistent & (isr.isr_min_abs >= 0.5)
    isr["isr_direction"] = np.where(isr.isr_dTE_mean > 0, "TE_up", "TE_down")
    ov = cores[cores.in_core].merge(isr, left_on="gene_id", right_index=True, how="left")
    ov["tet2_direction"] = np.where(ov.log2fc > 0, "up_in_TET2", "down_in_TET2")
    for c in D.columns:
        if cmeta.context_class[c] != "ISR engaged" and species[cmeta.study_id[c]] == "mouse":
            ov["dTE|" + c.split("|", 1)[1]] = ov.gene_id.map(D[c])
    ov.to_csv(out / "intersection_gene_overlap.csv", index=False)
    rows = []
    for core, d in cores.groupby("core"):
        uni = d.gene_id[d.gene_id.isin(isr.index[isr.isr_n_contexts == len(engaged)])]
        a, b = set(d.gene_id[d.in_core]) & set(uni), set(isr.index[isr.isr_core]) & set(uni)
        k = len(a & b)
        p = float(sps.hypergeom.sf(k - 1, len(uni), len(a), len(b))) if len(a) and len(b) else np.nan
        rows.append(dict(tet2_core=core, universe=len(uni), n_tet2_core=len(a), n_isr_core=len(b), overlap=k,
                         expected=len(a) * len(b) / max(len(uni), 1), p_hypergeom=p,
                         genes=";".join(sorted(sym.reindex(list(a & b)).dropna()))))
    counts = pd.DataFrame(rows)
    counts.to_csv(out / "intersection_overlap_counts.csv", index=False)
    return dict(pairs=pairs, sets=sets, cores=cores, overlap=ov, counts=counts, null_q95=q95)


def figure(pairs: pd.DataFrame, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    order = ["ISR engaged", "stress, ISR gate not passed", "ISR-blocked control", "low depth"]
    ctx = (pairs.drop_duplicates("isr_context").assign(o=lambda d: d.context_class.map(order.index))
           .sort_values(["o", "isr_context"]))
    M = pairs.pivot(index="isr_context", columns="tet2_unit", values="partial_spearman").loc[ctx.isr_context]
    sm = pairs.pivot(index="isr_context", columns="tet2_unit", values="species_matched").loc[ctx.isr_context]
    fig, ax = plt.subplots(figsize=(0.62 * M.shape[1] + 4.2, 0.36 * M.shape[0] + 2.6))
    im = ax.imshow(M.values, cmap="RdBu_r", vmin=-0.3, vmax=0.3, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6,
                        fontweight="bold" if sm.values[i, j] else "normal", color="#0b0b0b")
    ax.set_xticks(range(M.shape[1]), [c.replace("|", "\n") for c in M.columns], fontsize=6.5, rotation=90)
    ax.set_yticks(range(M.shape[0]), [f"{c.split('|', 1)[1]}  [{k}]" for c, k in zip(ctx.isr_context, ctx.context_class)], fontsize=6.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("partial Spearman", fontsize=7)
    ax.set_title("TET2 RNA change (columns) vs stress-induced TE change (rows), gene-wise\n"
                 "partial rank correlation controlling for expression, CDS length, GC, 3'UTR length; bold = same species",
                 fontsize=9, loc="left")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return path
