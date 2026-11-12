"""Effect-scale arithmetic: the mathematical core of the project.

Definitions (gene g, genotype k in {WT, TET2}, ISR condition s in {basal, high}):
    S = M * TE                              relative synthesis potential (arbitrary units)
    dTE[g,k] = log2(TE[g,k,high] / TE[g,k,basal])
    I_TE[g]  = dTE[g,TET2] - dTE[g,WT]
    I_S[g]   = dlog2S[g,TET2] - dlog2S[g,WT]

All inputs are linear-scale, positive arrays unless a name starts with ``log2``.
Nothing here is a measured protein abundance or an absolute flux.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def dte(te_basal, te_high):
    return np.log2(np.asarray(te_high, float) / np.asarray(te_basal, float))


def interaction_te(te_wt_basal, te_wt_high, te_tet2_basal, te_tet2_high):
    return dte(te_tet2_basal, te_tet2_high) - dte(te_wt_basal, te_wt_high)


def synthesis(m, te):
    return np.asarray(m, float) * np.asarray(te, float)


def interaction_log_s(s_wt_basal, s_wt_high, s_tet2_basal, s_tet2_high):
    """Log-scale difference-in-differences of synthesis potential."""
    f = lambda a: np.log2(np.asarray(a, float))
    return (f(s_tet2_high) - f(s_tet2_basal)) - (f(s_wt_high) - f(s_wt_basal))


def interaction_linear_s(s_wt_basal, s_wt_high, s_tet2_basal, s_tet2_high):
    """Linear-scale difference-in-differences (arbitrary common units)."""
    a = lambda x: np.asarray(x, float)
    return (a(s_tet2_high) - a(s_tet2_basal)) - (a(s_wt_high) - a(s_wt_basal))


def aggregate_log_change(m_basal, dte_log2, m_high=None):
    """log2 change of a *sum* of synthesis potentials over a gene set.

    log2( sum_g M_high[g] * 2^dTE[g] / sum_g M_basal[g] ), with basal TE set to 1
    on an arbitrary common scale.  Because the weights are the mRNA abundances,
    two genotypes with identical per-gene dTE but different M give different
    aggregate changes: a pure substrate/compositional effect.
    """
    m_basal = np.asarray(m_basal, float)
    m_high = m_basal if m_high is None else np.asarray(m_high, float)
    return np.log2(np.sum(m_high * 2.0 ** np.asarray(dte_log2, float)) / np.sum(m_basal))


def two_by_two_example() -> dict[str, pd.DataFrame]:
    """Worked 2x2 example (brief section L4).

    Gene A: TET2 loss doubles basal mRNA; ISR doubles TE identically in both genotypes.
    Gene B: TET2 loss leaves mRNA unchanged; ISR halves TE identically.
    Pathway {A,B}: sum of the two.
    """
    rows = []
    m = {"A": {"WT": 10.0, "TET2": 20.0}, "B": {"WT": 40.0, "TET2": 40.0}}
    te = {"A": {"basal": 1.0, "high": 2.0}, "B": {"basal": 1.0, "high": 0.5}}
    for g in m:
        for k in ("WT", "TET2"):
            for s in ("basal", "high"):
                rows.append(dict(gene=g, genotype=k, isr=s, M=m[g][k], TE=te[g][s],
                                 S=m[g][k] * te[g][s]))
    cells = pd.DataFrame(rows)
    path = cells.groupby(["genotype", "isr"], as_index=False)[["M", "S"]].sum()
    path.insert(0, "gene", "pathway{A,B}")
    path["TE"] = path["S"] / path["M"]  # abundance-weighted mean TE
    cells = pd.concat([cells, path], ignore_index=True)

    out = []
    for g, d in cells.groupby("gene", sort=False):
        v = lambda col, k, s: float(d[(d.genotype == k) & (d.isr == s)][col].iloc[0])
        out.append(dict(
            unit=g,
            dTE_WT=dte(v("TE", "WT", "basal"), v("TE", "WT", "high")),
            dTE_TET2=dte(v("TE", "TET2", "basal"), v("TE", "TET2", "high")),
            I_TE=interaction_te(v("TE", "WT", "basal"), v("TE", "WT", "high"),
                                v("TE", "TET2", "basal"), v("TE", "TET2", "high")),
            dS_linear_WT=v("S", "WT", "high") - v("S", "WT", "basal"),
            dS_linear_TET2=v("S", "TET2", "high") - v("S", "TET2", "basal"),
            I_S_linear=interaction_linear_s(v("S", "WT", "basal"), v("S", "WT", "high"),
                                            v("S", "TET2", "basal"), v("S", "TET2", "high")),
            I_S_log2=interaction_log_s(v("S", "WT", "basal"), v("S", "WT", "high"),
                                       v("S", "TET2", "basal"), v("S", "TET2", "high")),
        ))
    return {"cells": cells, "summary": pd.DataFrame(out).round(4)}
