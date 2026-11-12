"""Differential translation efficiency from paired RNA + Ribo counts (brief section C).

Model per context, on log2-CPM with assay-specific size factors:
    y ~ 1 + stress + ribo + stress:ribo
RNA change  = stress;  Ribo-occupancy change = stress + stress:ribo;  dTE = stress:ribo.
This is the same contrast structure as DESeq2 `deltaTE`.  Ribosome occupancy is not a
synthesis rate: elongation changes and uORF footprints are not separated here.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..stats import lm_contrasts, log_cpm, random_effects_meta, size_factors


def context_dte(counts: pd.DataFrame, samples: pd.DataFrame, min_mean: float = 10.0,
                group_col: str = "condition", ref: str = "control") -> pd.DataFrame:
    s = samples.set_index("sample_id").loc[counts.columns]
    rna, ribo = (s.assay == "rna").values, (s.assay == "ribo").values
    keep = (counts.loc[:, rna].mean(axis=1) >= min_mean) & (counts.loc[:, ribo].mean(axis=1) >= min_mean)
    c = counts[keep].values
    y = np.empty(c.shape)
    for m in (rna, ribo):  # normalise each assay separately
        y[:, m] = log_cpm(c[:, m], size_factors(c[:, m]))
    st = (s[group_col] != ref).values.astype(float)
    rb = ribo.astype(float)
    X = np.column_stack([np.ones(len(s)), st, rb, st * rb])
    res = lm_contrasts(y, X, {"rna": np.array([0, 1.0, 0, 0]), "ribo": np.array([0, 1.0, 0, 1.0]),
                              "dte": np.array([0, 0, 0, 1.0])})
    out = pd.DataFrame(dict(gene_id=counts.index[keep].values))
    for k, name in (("rna", "dRNA"), ("ribo", "dRibo"), ("dte", "dTE")):
        r = res[k]
        out[f"{name}_log2fc"], out[f"{name}_se"] = r.log2fc.values, r.se.values
        out[f"{name}_ci_lo"], out[f"{name}_ci_hi"] = r.ci_lo.values, r.ci_hi.values
        out[f"{name}_fdr"] = r.fdr.values
    out["n_bio_reps_per_arm"] = int(s.groupby(["assay", group_col]).bio_rep.nunique().min())
    return out


def library_qc(counts: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    """QC available from processed counts only. Read-level QC (footprint length, periodicity,
    P-site offsets, CDS fraction, rRNA) needs raw reads and is reported as not available."""
    rows = []
    for (cid, assay), s in samples.groupby(["context_id", "assay"]):
        c = counts[s.sample_id]
        lc = log_cpm(c.values)
        cor = np.corrcoef(lc.T)
        cond = s.condition.to_numpy(dtype=object)
        same = cond[:, None] == cond[None, :]
        off = ~np.eye(len(s), dtype=bool)
        rows.append(dict(context_id=cid, assay=assay, n_libraries=len(s),
                         median_depth=float(c.sum().median()),
                         genes_detected=int((c.mean(axis=1) >= 10).sum()),
                         min_within_condition_r=float(cor[same & off].min()),
                         footprint_length_qc="not available (processed counts only)",
                         periodicity_qc="not available (processed counts only)",
                         rrna_fraction="not available (processed counts only)"))
    return pd.DataFrame(rows)


def build_signatures(ctx: pd.DataFrame, samples: pd.DataFrame, counts: pd.DataFrame) -> pd.DataFrame:
    """Long table: one row per gene x context with RNA, Ribo and TE changes kept separate."""
    out = []
    for cid, s in samples.groupby("context_id"):
        d = context_dte(counts[s.sample_id], s)
        d.insert(1, "context_id", cid)
        out.append(d)
    sig = pd.concat(out, ignore_index=True).merge(ctx, on="context_id")
    sig["provenance"] = sig.study_id
    sig["evidence"] = "observed"
    sig["control_evidence"] = "matched unstressed control, same study and cell type"
    return sig


def pooled_signature(sig: pd.DataFrame, contexts=None, per_unit_intensity: bool = False) -> pd.DataFrame:
    """Random-effects pooled dTE per gene across contexts. `pred_sd` (includes between-context
    tau) is the uncertainty to carry when transferring the signature to a NEW context."""
    if contexts is not None:
        sig = sig[sig.context_id.isin(contexts)]
    scale = sig.intensity if per_unit_intensity else 1.0
    d = sig.assign(e=sig.dTE_log2fc / scale, s=sig.dTE_se / scale)
    est = d.pivot(index="gene_id", columns="context_id", values="e")
    se = d.pivot(index="gene_id", columns="context_id", values="s")
    m = random_effects_meta(est.values, se.values)
    m.insert(0, "gene_id", est.index.values)
    return m.rename(columns={"estimate": "dTE"})


def genotype_by_isr_dte(counts: pd.DataFrame, samples: pd.DataFrame, min_mean: float = 10.0) -> pd.DataFrame:
    """Direct I_TE from a full genotype x ISR x assay design (three-way interaction).

    This is the ONLY estimator in the project that identifies I_TE.  It needs paired RNA and
    Ribo libraries for all four genotype x condition cells with biological replicates.
    """
    s = samples.set_index("sample_id").loc[counts.columns]
    rna, ribo = (s.assay == "rna").values, (s.assay == "ribo").values
    keep = (counts.loc[:, rna].mean(axis=1) >= min_mean) & (counts.loc[:, ribo].mean(axis=1) >= min_mean)
    c = counts[keep].values
    y = np.empty(c.shape)
    for m in (rna, ribo):
        y[:, m] = log_cpm(c[:, m], size_factors(c[:, m]))
    g = (s.genotype == "TET2").values.astype(float)
    st = (s.condition != "control").values.astype(float)
    rb = ribo.astype(float)
    X = np.column_stack([np.ones(len(s)), g, st, rb, g * st, g * rb, st * rb, g * st * rb])
    e = np.eye(8)
    res = lm_contrasts(y, X, {"dTE_WT": e[6], "dTE_TET2": e[6] + e[7], "I_TE": e[7],
                              "I_RNA": e[4]})
    out = pd.DataFrame(dict(gene_id=counts.index[keep].values))
    for k, r in res.items():
        out[k], out[f"{k}_se"], out[f"{k}_fdr"] = r.log2fc.values, r.se.values, r.fdr.values
    return out
