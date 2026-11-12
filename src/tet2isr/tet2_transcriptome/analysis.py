"""Replicate-aware basal TET2 contrasts within matched compartments (brief section B)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..data_ingest.orthologs import map_counts
from ..stats import lm_contrasts, log_cpm, random_effects_meta, size_factors


def pseudobulk(counts: pd.DataFrame, cells: pd.DataFrame, min_cells: int = 10):
    """Sum single-cell counts per (sample, compartment). Cells are never treated as replicates."""
    keys = cells.set_index("cell_id").loc[counts.columns]
    grp = keys["sample_id"] + "|" + keys["compartment"]
    n = grp.value_counts()
    keep = n[n >= min_cells].index
    pb = counts.T.groupby(grp.values).sum().T[sorted(keep)]
    meta = pd.DataFrame(dict(pb_id=pb.columns))
    meta[["sample_id", "compartment"]] = meta.pb_id.str.split("|", expand=True)
    meta["n_cells"] = meta.pb_id.map(n).values
    return pb, meta


def composition(cells: pd.DataFrame, samples: pd.DataFrame) -> pd.DataFrame:
    frac = (cells.groupby("sample_id")["compartment"].value_counts(normalize=True)
            .rename("fraction").reset_index())
    return frac.merge(samples[["sample_id", "genotype"]], on="sample_id")


def genotype_de(counts: pd.DataFrame, genotype: pd.Series, min_mean: float = 5.0) -> pd.DataFrame:
    """log2(TET2/WT) per gene. One column of `counts` = one independent biological replicate."""
    g = (genotype.values == "TET2").astype(float)
    if min(g.sum(), (1 - g).sum()) < 2:
        raise ValueError("need >=2 independent biological replicates per genotype")
    keep = counts.mean(axis=1) >= min_mean
    c = counts[keep]
    y = log_cpm(c.values, size_factors(c.values))
    X = np.column_stack([np.ones(len(g)), g])
    res = lm_contrasts(y, X, {"tet2_vs_wt": np.array([0.0, 1.0])})["tet2_vs_wt"]
    res.insert(0, "gene_id", c.index.values)
    res["mean_log2cpm_WT"] = y[:, g == 0].mean(axis=1)
    res["mean_log2cpm_TET2"] = y[:, g == 1].mean(axis=1)
    res["n_WT"], res["n_TET2"] = int((1 - g).sum()), int(g.sum())
    return res


def run_all(studies: dict, target_compartment: str = "HSC") -> dict:
    """Per-study, per-compartment contrasts + random-effects meta-analysis in the target compartment."""
    per, qc, comp, expr = [], [], [], {}
    for sid, s in studies.items():
        counts, smp = s["counts"], s["samples"]
        info = {}
        if s["orthologs"] is not None:
            counts, info = map_counts(counts, s["orthologs"])
        if s["cells"] is not None:
            comp.append(composition(s["cells"], smp).assign(study_id=sid))
            pb, meta = pseudobulk(counts, s["cells"])
            meta = meta.merge(smp[["sample_id", "genotype", "donor"]], on="sample_id")
            units = [(c, pb[m.pb_id], m.genotype.reset_index(drop=True), m.n_cells.values)
                     for c, m in meta.groupby("compartment")]
        else:
            units = [(c, counts[m.sample_id], m.genotype.reset_index(drop=True), None)
                     for c, m in smp.groupby("compartment")]
        for compartment, cnt, geno, ncell in units:
            de = genotype_de(cnt, geno).assign(study_id=sid, compartment=compartment)
            per.append(de)
            lc = log_cpm(cnt.values)
            cor = np.corrcoef(lc.T)
            qc.append(dict(study_id=sid, compartment=compartment, n_samples=cnt.shape[1],
                           n_WT=int((geno == "WT").sum()), n_TET2=int((geno == "TET2").sum()),
                           median_depth=float(cnt.sum().median()),
                           min_replicate_r=float(cor.min()), genes_tested=len(de),
                           median_cells_per_pseudobulk=(float(np.median(ncell)) if ncell is not None else np.nan),
                           **info))
            if compartment == target_compartment:
                expr[sid] = dict(counts=cnt, genotype=geno)
    per = pd.concat(per, ignore_index=True)
    tgt = per[per.compartment == target_compartment]
    est = tgt.pivot(index="gene_id", columns="study_id", values="log2fc")
    se = tgt.pivot(index="gene_id", columns="study_id", values="se")
    meta = random_effects_meta(est.values, se.values)
    meta.insert(0, "gene_id", est.index.values)
    meta = meta.rename(columns={"estimate": "log2fc"})
    meta["n_studies_same_sign"] = (np.sign(est.values) == np.sign(meta.log2fc.values)[:, None]).sum(axis=1)
    meta["tau_note"] = np.where(meta.k < 3, "tau poorly identified (k<3 studies)", "")
    return dict(per_study=per, meta=meta, qc=pd.DataFrame(qc), expr=expr,
                composition=(pd.concat(comp) if comp else pd.DataFrame()))
