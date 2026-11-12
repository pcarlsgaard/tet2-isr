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
    lab = ["basal TET2\ncontrast", "ISR\ntranscription", "ISR translation\n(RNA+Ribo)", "direct TET2\nx ISR TE"]
    m = manifest[manifest.source_class.isin(["tet2_hspc", "isr_translation", "isr_transcription"])].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.6, 0.36 * len(m) + 1.3))
    for i, r in m.iterrows():
        for j, c in enumerate(caps):
            v = str(r[c])
            colr = WT if v == "yes" else "#ecebe7" if v.startswith("no") and "(" not in v or v == "no" else "#f3d9a4"
            ax.add_patch(plt.Rectangle((j + 0.04, i + 0.08), 0.92, 0.84, color=colr))
            ax.text(j + 0.5, i + 0.5, v.replace("claimed_unverified", "claimed, unverified").replace(" (", "\n("),
                    ha="center", va="center", fontsize=6.5, color="white" if v == "yes" else INK)
        used = str(r.get("use_in_run", ""))
        ax.text(4.25, i + 0.5, used[:58], va="center", fontsize=6.5,
                color=MUTED if used.lower().startswith(("not", "excluded")) else INK)
    ax.set_xlim(0, 6.9); ax.set_ylim(len(m), 0)
    ax.set_xticks(list(np.arange(4) + 0.5) + [4.25], lab + ["\nuse in this run"])
    ax.get_xticklabels()[-1].set_ha("left")
    ax.xaxis.tick_top(); ax.tick_params(length=0)
    ax.set_yticks(np.arange(len(m)) + 0.5, [f"{str(a)[:22]}  ({s})" for a, s in zip(m.accession, m.species)])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Evidence map: what each public dataset can support (amber = limited or qualified; see manifest)", pad=34)
    return _save(fig, path)


def fig2_basal(m0: pd.DataFrame, path):
    progs = list(m0.groupby("program").log2fc.mean().sort_values().index)
    studies = list(m0.study_id.unique())
    fig, ax = plt.subplots(figsize=(6.6, (0.2 + 0.11 * len(studies)) * len(progs) + 1.4))
    for k, (s, c) in enumerate(zip(studies, (WT, KO, ACC, "#eda100"))):
        d = m0[m0.study_id == s].set_index("program").reindex(progs)
        y = np.arange(len(progs)) + (k - (len(studies) - 1) / 2) * (0.7 / len(studies))
        ax.errorbar(d.log2fc, y, xerr=[d.log2fc - d.ci_lo, d.ci_hi - d.log2fc], fmt="o", color=c, ms=4,
                    lw=1.3, label=f"{s} (n={int(d.n_WT.dropna().iloc[0])} vs {int(d.n_TET2.dropna().iloc[0])})")
    ax.axvline(0, color=GREY, lw=1)
    ax.set_yticks(range(len(progs)), [p.replace("HALLMARK_", "H: ").lower() if p.startswith("HALLMARK") else p for p in progs])
    ax.set_xlabel("program score, log2(TET2 / WT), 95% CI over biological replicates")
    ax.set_title("M0: basal RNA programs by genotype (target population of each study)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, fontsize=7)
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
    """Every context predicted by the pooled-signature model trained without that context's study."""
    from ..modeling.models import make_model
    cids = list(ctx.index)
    ncol = min(4, len(cids))
    nrow = int(np.ceil(len(cids) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.4 * ncol, 3.3 * nrow), squeeze=False)
    for ax in axes.ravel()[len(cids):]:
        ax.axis("off")
    for ax, c in zip(axes.ravel(), cids):
        tr = ctx[ctx.study_id != ctx.loc[c, "study_id"]]
        m = make_model("M1b_signature+dose", genes, state_cols).fit(sig[sig.context_id.isin(tr.index)], genes, tr)
        o = sig[sig.context_id == c].set_index("gene_id").dTE_log2fc
        p = m.predict(genes, ctx.loc[[c]])[c].reindex(o.index)
        o, p = o - o.median(), p - p.median()
        ax.scatter(p, o, s=3, color=WT, alpha=0.35, rasterized=True)
        ax.axhline(0, color=GREY, lw=0.6); ax.axvline(0, color=GREY, lw=0.6)
        ax.set_title(f"{c.split('|', 1)[1][:34]}\nr = {np.corrcoef(o, p)[0, 1]:.2f}", fontsize=8)
        ax.set_xlabel("predicted dTE", fontsize=8); ax.set_ylabel("observed dTE", fontsize=8)
    fig.suptitle("Held-out study: pooled dTE signature from the other studies vs observed (dose-blind)",
                 x=0.01, ha="left", fontsize=10)
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


def fig6_genotype(m1p, hp, path, m2_note=None):
    d = m1p[m1p.isr_intensity == "moderate"]
    progs = sorted(d.program.unique())
    studies = list(d.tet2_study.unique())
    lab = [p.replace("HALLMARK_", "H: ").lower() if p.startswith("HALLMARK") else p for p in progs]
    two = hp is not None and len(hp)
    fig, axs = plt.subplots(1, 2 if two else 1, figsize=(10.5 if two else 6.6, 0.3 * len(progs) + 1.8), sharey=True, squeeze=False)
    a = axs[0, 0]
    for k, (s, c) in enumerate(zip(studies, (WT, KO, ACC, "#eda100"))):
        x = d[d.tet2_study == s].set_index("program").reindex(progs)
        y = np.arange(len(progs)) + (k - (len(studies) - 1) / 2) * (0.7 / len(studies))
        a.errorbar(x.I_S_log2_program, y, xerr=[x.I_S_log2_program - x.I_S_log2_program_lo,
                                                x.I_S_log2_program_hi - x.I_S_log2_program],
                   fmt="o", ms=4, lw=1.2, color=c, label=s)
    a.axvline(0, color=GREY, lw=1); a.set_yticks(range(len(progs)), lab)
    a.set_xlabel("program-level I_S, log2 (95% MC interval)")
    a.set_title("M1: substrate-only interaction per TET2 study\n(gene-level I_TE = 0 by assumption)")
    a.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4)
    if m2_note:
        a.text(0, -0.3, m2_note, transform=a.transAxes, fontsize=7, color=MUTED, va="top", wrap=True)
    if two:
        b = axs[0, 1]
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
        b.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    return _save(fig, path)


def fig8_reference(refchk: pd.DataFrame, eng: pd.DataFrame, path):
    """Reference uORF-target dTE in every context, grouped by role: the positive control."""
    r = refchk.copy()
    r["label"] = r.context_id.str.split("|", n=1).str[1]
    e = eng.set_index("context_id").engaged if eng is not None and "engaged" in eng else pd.Series(dtype=bool)
    grp = lambda row: ("ISR blocked (ISRIB / eIF2a-S51A / Hri-KO)" if row.role == "specificity_control" else
                       "excluded: low depth" if row.role == "excluded_low_depth" else
                       "engaged (defines signature)" if e.get(row.context_id, False) else "stress, gate not passed")
    r["group"] = r.apply(grp, axis=1)
    order = ["engaged (defines signature)", "stress, gate not passed", "ISR blocked (ISRIB / eIF2a-S51A / Hri-KO)",
             "excluded: low depth"]
    col = dict(zip(order, (WT, "#eda100", KO, GREY)))
    labs = (r.groupby(["group", "label"]).dTE_minus_context_median.mean().reset_index()
            .assign(o=lambda d: d.group.map(order.index)).sort_values(["o", "dTE_minus_context_median"], ascending=[True, False]))
    fig, ax = plt.subplots(figsize=(7.6, 0.33 * len(labs) + 1.6))
    seen = set()
    for i, row in enumerate(labs.itertuples()):
        x = r[r.label == row.label].dTE_minus_context_median
        ax.scatter(x, np.full(len(x), i), s=14, color=col[row.group], alpha=0.75,
                   label=row.group if row.group not in seen else None)
        ax.plot([x.mean()] * 2, [i - 0.3, i + 0.3], color=INK, lw=1.5)
        seen.add(row.group)
    ax.axvline(0, color=GREY, lw=1)
    ax.set_yticks(range(len(labs)), labs.label, fontsize=7); ax.invert_yaxis()
    ax.set_xlabel("dTE of reference uORF targets, log2 (relative to context median); bar = set mean")
    ax.set_title("Positive control: do reference ISR targets gain TE in each context?")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
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


def gene_table_public(out: Path, genes: pd.DataFrame, programs: pd.DataFrame, ctx: pd.DataFrame) -> pd.DataFrame:
    """One row per gene: observed basal TET2 contrast, transferred dTE with between-context spread,
    and the M1 substrate-only linear interaction. No I_TE estimate is made."""
    sym = genes.set_index("gene_id").symbol
    meta = pd.read_parquet(out / "tet2_gene_contrasts_meta.parquet").set_index("gene_id")
    sig = pd.read_parquet(out / "isr_translation_signatures.parquet")
    sig = sig[sig.context_id.isin(ctx.index)]
    Y = sig.pivot(index="gene_id", columns="context_id", values="dTE_log2fc")
    Y = Y - Y.median()
    m1 = pd.read_parquet(out / "m1_gene_predictions.parquet")
    m1 = m1[m1.isr_intensity == "moderate"]
    nz = (m1.I_S_linear_lo > 0) | (m1.I_S_linear_hi < 0)
    g1 = m1.assign(pos=nz & (m1.I_S_linear > 0), neg=nz & (m1.I_S_linear < 0)).groupby("gene_id").agg(
        dTE_transferred=("dTE_transferred", "first"), dTE_pred_sd=("dTE_pred_sd", "first"),
        I_S_linear_median=("I_S_linear", "median"), n_tet2_studies=("tet2_study", "nunique"),
        n_pos=("pos", "sum"), n_neg=("neg", "sum"))
    t = pd.DataFrame(dict(symbol=sym)).join(g1, how="inner")
    t["n_isr_contexts"] = Y.notna().sum(axis=1).reindex(t.index)
    t["n_contexts_dTE_up"] = (Y > 0).sum(axis=1).reindex(t.index)
    t["dTE_min"], t["dTE_max"] = Y.min(axis=1).reindex(t.index), Y.max(axis=1).reindex(t.index)
    t["dTE_lo"], t["dTE_hi"] = t.dTE_transferred - 1.96 * t.dTE_pred_sd, t.dTE_transferred + 1.96 * t.dTE_pred_sd
    t["dTE_consistent"] = (t.n_isr_contexts == len(ctx)) & ((t.n_contexts_dTE_up == len(ctx)) | (t.n_contexts_dTE_up == 0)) \
        & (np.minimum(Y.abs().min(axis=1).reindex(t.index), np.inf) >= 0.5)
    for c in ("log2fc", "se", "fdr", "k", "n_studies_same_sign", "replicated"):
        t[f"basal_tet2_{c}"] = meta[c].reindex(t.index)
    t["substrate_interaction_supported"] = (np.maximum(t.n_pos, t.n_neg) >= 3) & (np.minimum(t.n_pos, t.n_neg) == 0)
    t["programs"] = programs.groupby("gene_id").program.agg(lambda x: ";".join(sorted(x))).reindex(t.index).fillna("")
    t.index.name = "gene_id"
    t.reset_index().to_csv(out / "gene_predictions.csv", index=False)
    return t


def priority_table_public_run(out: Path, genes, programs, ctx, prior_path="configs/prior_hypotheses.tsv",
                              dest="results/prediction_priority_table.csv") -> pd.DataFrame:
    """Final priority table: prior hypotheses joined with what the public data do and do not show."""
    t = gene_table_public(out, genes, programs, ctx)
    m1p = pd.read_csv(out / "m1_program_predictions.csv")
    m1p = m1p[m1p.isr_intensity == "moderate"]
    m0 = pd.read_csv(out / "m0_program_contrasts.csv")
    nctx, src_isr = len(ctx), ";".join(sorted(ctx.study_id.unique()))
    src_tet2 = ";".join(sorted(m1p.tet2_study.unique()))
    direction = lambda v, sd: "up" if v - 1.96 * sd > 0 else "down" if v + 1.96 * sd < 0 else "not distinguishable from 0"
    rows = []

    def gene_row(sym, mech, alt, readout, note=""):
        r = t[t.symbol == sym]
        if r.empty:
            return
        r = r.iloc[0]
        d = direction(r.dTE_transferred, r.dTE_pred_sd)
        rows.append(dict(
            gene_or_program=sym, level="gene",
            predicted_direction_WT=d, predicted_direction_TET2=d + " (same by M1 assumption)",
            dTE_transferred_log2=round(r.dTE_transferred, 2),
            dTE_prediction_interval=f"[{r.dTE_lo:.2f}, {r.dTE_hi:.2f}]",
            dTE_range_across_contexts=f"[{r.dTE_min:.2f}, {r.dTE_max:.2f}] in {int(r.n_isr_contexts)}/{nctx} contexts",
            basal_TET2_log2fc_meta=round(r.basal_tet2_log2fc, 2) if pd.notna(r.basal_tet2_log2fc) else np.nan,
            basal_TET2_evidence=("replicated" if r.basal_tet2_replicated is True or r.basal_tet2_replicated == 1 else
                                 f"meta FDR {r.basal_tet2_fdr:.2f}; same sign in "
                                 f"{int(r.basal_tet2_n_studies_same_sign)}/{int(r.basal_tet2_k)} studies; "
                                 + ("fails leave-heaviest-study-out" if r.basal_tet2_fdr < 0.05 else "not significant")
                                 if pd.notna(r.basal_tet2_fdr) else "not measured"),
            predicted_I_TE="not identifiable from public data (M1 assumes 0; M2 not supported)",
            predicted_change_synthesis_potential=(
                f"substrate-only linear interaction {'supported' if r.substrate_interaction_supported else 'not supported'}"
                f" ({int(max(r.n_pos, r.n_neg))}/{int(r.n_tet2_studies)} TET2 studies with interval excluding 0)"),
            evidence_category="transferred (dTE); observed (basal RNA)", source_ids=f"ISR: {src_isr} | TET2: {src_tet2}",
            stress_dose_time="moderate ISR, early translational phase; contexts span 1-16 h",
            model_stability=("consistent sign, |dTE|>=0.5 in all engaged contexts" if r.dTE_consistent else
                             f"dTE up in {int(r.n_contexts_dTE_up)}/{int(r.n_isr_contexts)} contexts") + note,
            plausible_mechanism=mech, alternative_explanation=alt, recommended_discriminating_readout=readout))

    prior = pd.read_csv(prior_path, sep="\t")
    pr = lambda key: prior[prior.gene_or_program.str.startswith(key)].iloc[0]
    circ = "; member of the reference set used to gate contexts (partly circular)"
    for sym, key in (("ATF4", "ATF4 "), ("ATF5", "ATF5"), ("DDIT3", "DDIT3"), ("PPP1R15A", "PPP1R15A")):
        q = pr(key)
        gene_row(sym, q.plausible_mechanism, q.alternative_explanation, q.recommended_discriminating_readout, circ)
    for sym in ("IFRD1", "SLC35A4", "PTP4A1", "PPP1R15B", "UCP2"):
        gene_row(sym, "uORF-dependent translation resistant to low ternary complex", "5'UTR isoform differences between cell types",
                 "5'UTR reporter; polysome qPCR", circ)
    new = t[t.dTE_consistent & ~t.programs.str.contains("isr_uorf_targets")].assign(a=lambda d: d.dTE_transferred.abs()).nlargest(12, "a")
    for sym in new.symbol:
        gene_row(sym, "unknown; not in the reference set", "three contexts only: chance agreement likely for some genes",
                 "replicate in an independent paired RNA/Ribo ISR dataset before any follow-up",
                 "; data-derived candidate, not a prior ISR target")
    for sym in ("MCL1", "BCL2", "BCL2L1", "PMAIP1", "BCL2L11", "BBC3"):
        q = pr("Apoptotic priming")
        gene_row(sym, q.plausible_mechanism, q.alternative_explanation, q.recommended_discriminating_readout)
    for sym in ("EIF2AK1", "TET2", "NCOA4"):
        gene_row(sym, "context gene for the HRI / iron hypotheses", "none tested", "see heme/iron row in configs/prior_hypotheses.tsv")
    for p_, d in m1p.groupby("program"):
        nz = ((d.I_S_log2_program_lo > 0) | (d.I_S_log2_program_hi < 0))
        b = m0[m0.program == p_]
        rows.append(dict(
            gene_or_program=p_, level="program",
            predicted_direction_WT=("up" if d.dlog2S_WT.mean() > 0.2 else "down" if d.dlog2S_WT.mean() < -0.2 else "flat"),
            predicted_direction_TET2=("up" if d.dlog2S_TET2.mean() > 0.2 else "down" if d.dlog2S_TET2.mean() < -0.2 else "flat"),
            dTE_transferred_log2=round(d.dlog2S_WT.mean(), 2),
            dTE_prediction_interval="program aggregate of synthesis potential; see m1_program_predictions.csv",
            dTE_range_across_contexts="",
            basal_TET2_log2fc_meta=round(b.log2fc.mean(), 2),
            basal_TET2_evidence=f"program score p<0.05 in {int((b.pvalue < 0.05).sum())}/{len(b)} TET2 studies; "
                                f"signs +{int((b.log2fc > 0).sum())}/-{int((b.log2fc < 0).sum())}",
            predicted_I_TE="not identifiable from public data",
            predicted_change_synthesis_potential=(
                f"program I_S (log2) range [{d.I_S_log2_program.min():.3f}, {d.I_S_log2_program.max():.3f}] across TET2 studies; "
                f"interval excludes 0 in {int(nz.sum())}/{len(d)}"),
            evidence_category="transferred", source_ids=f"ISR: {src_isr} | TET2: {src_tet2}",
            stress_dose_time="moderate ISR, early translational phase",
            model_stability="substrate-only; sign " + ("consistent" if abs(np.sign(d.I_S_log2_program).sum()) == len(d) else "inconsistent") + " across TET2 studies",
            plausible_mechanism="abundance-weighted aggregation of a common dTE (compositional)",
            alternative_explanation="cell-composition or batch differences between WT and TET2 libraries",
            recommended_discriminating_readout="set-level I_TE from a 2x2 paired RNA/Ribo design at matched p-eIF2a"))
    hyp = prior[prior.gene_or_program.str.startswith(("H0", "H1", "H2", "H3", "Heme", "Mitochondrial", "Basal ISR"))]
    for _, q in hyp.iterrows():
        rows.append(dict(gene_or_program=q.gene_or_program, level="hypothesis", predicted_direction_WT=q.predicted_direction_WT,
                         predicted_direction_TET2=q.predicted_direction_TET2, predicted_I_TE=q.predicted_I_TE,
                         predicted_change_synthesis_potential="not tested by public data",
                         evidence_category=q.evidence_category, source_ids=q.source_ids, stress_dose_time=q.stress_dose_time,
                         model_stability="no model fit", plausible_mechanism=q.plausible_mechanism,
                         alternative_explanation=q.alternative_explanation,
                         recommended_discriminating_readout=q.recommended_discriminating_readout))
    out_t = pd.DataFrame(rows)
    Path(dest).parent.mkdir(exist_ok=True)
    out_t.to_csv(dest, index=False)
    out_t.to_csv(out / "prediction_priority_table.csv", index=False)
    return out_t


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
    def rd(name, **k):
        try:
            return pd.read_csv(out / name, **k)
        except (FileNotFoundError, pd.errors.EmptyDataError):
            return None
    programs = pd.read_csv(root / "programs.tsv", sep="\t")
    programs.to_csv(out / "_programs.tsv", sep="\t", index=False)
    genes = pd.read_parquet(root / "genes.parquet")
    manifest = pd.read_csv("data/dataset_manifest.tsv", sep="\t")
    ctx = pd.read_csv(out / "isr_context_table.tsv", sep="\t").set_index("context_id")
    sig_all = pd.read_parquet(out / "isr_translation_signatures.parquet")
    sig = sig_all[sig_all.context_id.isin(ctx.index)]
    qc = pd.read_csv(out / "isr_library_qc.tsv", sep="\t")
    m0 = pd.read_csv(out / "m0_program_contrasts.csv")
    m1p = pd.read_csv(out / "m1_program_predictions.csv")
    ev = pd.read_csv(out / "heldout_metrics.csv")
    hp, hv = rd("holdout_program_level.csv"), rd("holdout_validation.csv")
    eng, refchk, spec = rd("isr_engagement.tsv", sep="\t"), rd("reference_set_check.csv"), rd("isr_specificity.csv")
    uneng, conc = rd("unengaged_context_check.csv"), rd("tet2_cross_study_spearman.csv")
    m2_missing = out / "m2_genotype_predictions.NOT_SUPPORTED.txt"
    m2_note = m2_missing.read_text().strip() if m2_missing.exists() else None
    ref = set(programs.loc[programs.program == cfg["reference_program"], "gene_id"])
    f = dict(evidence=fig1_evidence_map(manifest, figs / "fig1_evidence_map.png"),
             basal=fig2_basal(m0, figs / "fig2_basal_programs.png"),
             rna_te=fig3_rna_vs_te(sig_all, qc, ref, figs / "fig3_rna_vs_te_qc.png"),
             heldout=fig4_obs_pred(sig, ctx, genes, cfg["state_programs"], figs / "fig4_heldout_obs_vs_pred.png"),
             ablation=fig5_ablation(ev, figs / "fig5_ablation.png"),
             genotype=fig6_genotype(m1p, hp, figs / "fig6_genotype_by_isr.png"),
             hyp=fig7_hypotheses(figs / "fig7_hypotheses.png"))
    if refchk is not None:
        f["reference"] = fig8_reference(refchk, eng, figs / "fig8_reference_set.png")
    if synthetic:
        pri = priority_table_run(out, synthetic)
        pub = pd.read_csv("results/prediction_priority_table.csv") if Path("results/prediction_priority_table.csv").exists() \
            else priority_table_public()
    else:
        pub = pri = priority_table_public_run(out, genes, programs, ctx)
    metric_cols = ["pearson", "spearman", "rmse", "sign_concordance_sig", "top50_precision", "program_recovery_r",
                   "coverage_50", "coverage_80", "coverage_95"]
    summ = ev.groupby(["split", "model"])[metric_cols].mean().reset_index()
    summ.to_csv(out / "heldout_summary.csv", index=False)
    seen = (ev[(ev.split == "held_out_stressor") & ev.model.isin(["M1b_signature+dose", "M2_+cell_state"])]
            .groupby(["model", "stressor_class", "cell_type_in_training"])[["pearson", "rmse"]].mean().reset_index())
    banner = ("SYNTHETIC DEMONSTRATION. Every number, figure and gene ID (SYN...) below comes from simulated data "
              "with planted effects. It shows that the software works and what the design can and cannot "
              "identify. It says nothing about real TET2 biology." if synthetic else
              "PUBLIC-DATA RUN. Basal TET2 contrasts and ISR translation signatures are observed in the cited GEO "
              "series. Everything that combines them is a transfer-based prediction. No public dataset measures "
              "translation in TET2-deficient cells under ISR, so the genotype x ISR translation interaction (I_TE) "
              "is not estimated.")
    title = f"TET2 x ISR translatomic transfer: research report ({cfg['run_name']})"
    h = [f"<html><head><meta charset='utf-8'><title>{title}</title><style>{CSS}</style></head><body>",
         f"<h1>{title}</h1>", f"<p class='banner'>{banner}</p>",
         f"<p class='note'>Regenerate with <code>make {'smoke-test' if synthetic else 'public'}</code>. Methods: "
         "docs/methods.md. Limits: docs/assumptions_and_limits.md. Inputs and checksums: provenance.json in the run "
         "directory.</p>",
         "<h2>1. Datasets and evidence map</h2>", _img(f["evidence"]),
         "<p class='note'>Full fields, including how each was verified: data/dataset_manifest.tsv.</p>",
         _tbl(manifest[manifest.source_class.isin(["tet2_hspc", "isr_translation", "isr_transcription"])]
              [["accession", "species", "cell_type", "n_independent_bio_reps", "assay", "use_in_run", "limitations"]]),
         "<h2>2. Identifiability: the 2x2 example</h2>",
         "<p>Gene A doubles in mRNA with TET2 loss; ISR doubles its TE in both genotypes. Linear output differs "
         "(I_S_linear = 10 units) while I_TE = 0 and gene-level log I_S = 0. The two-gene pathway shows a non-zero "
         "<em>aggregate</em> log interaction purely because the abundance weights differ.</p>",
         _tbl(pd.read_csv(out / "example_2x2_summary.csv"), r=4),
         "<h2>3. Basal genotype contrasts (M0)</h2>", _img(f["basal"]),
         _tbl(pd.read_csv(out / "tet2_sample_qc.tsv", sep="\t").dropna(axis=1, how="all"))]
    if conc is not None:
        meta = pd.read_parquet(out / "tet2_gene_contrasts_meta.parquet")
        h += ["<p class='note'>Spearman correlation of gene-level log2(TET2/WT) between studies (target populations):</p>",
              _tbl(conc, r=2),
              f"<p>{int((meta.fdr < 0.05).sum())} genes reach meta-analysis FDR &lt; 0.05, but for "
              f"{int((meta[meta.fdr < 0.05].max_weight_share > 0.8).sum())} of them one study carries more than 80% of "
              f"the weight. {int(meta.replicated.sum())} genes are <em>replicated</em>: same sign in at least 3 studies and "
              "still significant with the heaviest study removed.</p>"]
    h += ["<h2>4. ISR translation signatures and QC</h2>"]
    if "reference" in f:
        h += [_img(f["reference"]),
              "<p class='note'>Engagement gate (reference-set mean dTE &ge; threshold and one-sided p &lt; alpha, "
              "configs/public.yml). Only engaged contexts define the transferred signature.</p>",
              _tbl(eng[["context_id", "n_ref_genes", "ref_mean_dTE", "ref_frac_positive", "ref_p_one_sided", "engaged"]])]
    if spec is not None and len(spec):
        h += ["<p class='note'>Each stress signature against its ISR-blocked counterpart:</p>", _tbl(spec)]
    h += [_img(f["rna_te"]),
          "<h2>5. Held-out prediction of dTE</h2>", _img(f["heldout"]), _img(f["ablation"]),
          "<p class='note'>Mean over held-out contexts. coverage_XX = fraction of held-out genes inside the "
          "nominal XX% predictive interval.</p>", _tbl(summ)]
    if len(seen):
        h += ["<p class='note'>Held-out stressor, by stressor class and whether the test cell type was seen:</p>", _tbl(seen)]
    if uneng is not None and len(uneng):
        h += ["<p class='note'>Negative control: agreement of the engaged-context signature with stress contexts that "
              "did not pass the gate:</p>", _tbl(uneng)]
    h += ["<p class='note'>Gene-family leakage check (sequence-only model, unseen genes):</p>",
          _tbl(pd.read_csv(out / "gene_family_leakage.csv")),
          "<p class='note'>M3 factor stability:</p>", _tbl(pd.read_csv(out / "m3_factor_stability.csv")),
          "<h2>6. Genotype x ISR predictions</h2>", _img(f["genotype"])]
    if m2_note:
        h.append(f"<p class='banner'>{m2_note}</p>")
    if hv is not None:
        h += ["<p class='note'>Gene-level comparison with the held-out direct 2x2 translatome (moderate dose):</p>",
              _tbl(hv[hv.isr_intensity == "moderate"][["intensity_source", "model", "tet2_study", "target",
                                                        "pearson", "rmse", "n_sig_obs", "top50_precision",
                                                        "coverage_95"]]),
              "<p class='note'>Program level:</p>", _tbl(hp)]
    else:
        h.append("<p class='banner'>No direct genotype x ISR translatome is available: nothing here validates a "
                 "genotype-specific translation effect.</p>")
    if synthetic:
        h += ["<h3>Run-specific priority table (top rows)</h3>", _tbl(pri, 16)]
    h += ["<h2>7. Prediction priority table, competing hypotheses and discriminating assays</h2>", _img(f["hyp"]),
          "<p class='note'>Full table: results/prediction_priority_table.csv</p>",
          _tbl(pub[[c for c in ("gene_or_program", "level", "predicted_direction_WT", "dTE_transferred_log2",
                                "dTE_prediction_interval", "dTE_range_across_contexts", "basal_TET2_log2fc_meta",
                                "basal_TET2_evidence", "predicted_I_TE", "predicted_change_synthesis_potential",
                                "evidence_category", "model_stability", "recommended_discriminating_readout")
                    if c in pub.columns]]),
          "</body></html>"]
    Path("reports").mkdir(exist_ok=True)
    p = Path("reports") / ("research_report_synthetic_demo.html" if synthetic else "research_report.html")
    p.write_text("\n".join(h))
    return p
