"""Shared data contract. Every ingest path (synthetic, GEO, in-house) must emit these tables.

<root>/genes.parquet            gene_id, family_id, seq_* (sequence features), epi_* (optional)
<root>/programs.tsv             program, gene_id
<root>/isr/contexts.tsv         CONTEXT_COLS
<root>/isr/samples.tsv          ISR_SAMPLE_COLS
<root>/isr/counts.parquet       index gene_id, columns sample_id (raw integer counts)
<root>/tet2/<study>/samples.tsv TET2_SAMPLE_COLS
<root>/tet2/<study>/counts.parquet   genes x samples (bulk) or genes x cells (single cell)
<root>/tet2/<study>/cells.tsv   (single cell only) cell_id, sample_id, compartment
<root>/tet2_isr_holdout/        optional direct genotype x ISR translatome; validation ONLY
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

CONTEXT_COLS = ["context_id", "study_id", "species", "cell_type", "stressor", "stressor_class",
                "kinase", "intensity", "duration_h"]
ISR_SAMPLE_COLS = ["sample_id", "context_id", "study_id", "condition", "assay", "bio_rep"]
TET2_SAMPLE_COLS = ["sample_id", "study_id", "species", "genotype", "allele_type", "compartment",
                    "donor"]
HOLDOUT_SAMPLE_COLS = ["sample_id", "genotype", "condition", "assay", "bio_rep"]
EVIDENCE = ("observed", "transferred", "speculative")


class ContractError(ValueError):
    pass


def require_cols(df: pd.DataFrame, cols, name: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ContractError(f"{name}: missing columns {missing}")


def check_counts(counts: pd.DataFrame, samples: pd.DataFrame, name: str, id_col="sample_id"):
    """Counts columns must match the sample sheet exactly, be non-negative and unique."""
    if counts.index.has_duplicates:
        raise ContractError(f"{name}: duplicated gene IDs")
    ids = list(samples[id_col])
    if len(set(ids)) != len(ids):
        raise ContractError(f"{name}: duplicated sample IDs")
    if set(ids) != set(counts.columns):
        raise ContractError(f"{name}: count columns and sample sheet differ: "
                            f"{sorted(set(ids) ^ set(counts.columns))[:5]}")
    if (counts.values < 0).any():
        raise ContractError(f"{name}: negative counts")
    return counts[ids]


def check_paired_rna_ribo(samples: pd.DataFrame) -> None:
    """Every (context, condition) needs >=2 biological replicates in BOTH assays."""
    require_cols(samples, ISR_SAMPLE_COLS, "isr samples")
    n = samples.groupby(["context_id", "condition", "assay"])["bio_rep"].nunique().unstack("assay")
    for a in ("rna", "ribo"):
        if a not in n or n[a].isna().any() or (n[a] < 2).any():
            raise ContractError(f"isr samples: assay '{a}' lacks >=2 biological replicates in "
                                f"some context/condition:\n{n}")


def load_isr(root: Path):
    root = Path(root)
    ctx = pd.read_csv(root / "isr/contexts.tsv", sep="\t")
    smp = pd.read_csv(root / "isr/samples.tsv", sep="\t")
    require_cols(ctx, CONTEXT_COLS, "contexts")
    check_paired_rna_ribo(smp)
    cnt = check_counts(pd.read_parquet(root / "isr/counts.parquet"), smp, "isr counts")
    return ctx, smp, cnt


def load_tet2_studies(root: Path) -> dict:
    out = {}
    for d in sorted((Path(root) / "tet2").iterdir()):
        smp = pd.read_csv(d / "samples.tsv", sep="\t")
        require_cols(smp, TET2_SAMPLE_COLS, f"tet2/{d.name} samples")
        cnt = pd.read_parquet(d / "counts.parquet")
        cells = pd.read_csv(d / "cells.tsv", sep="\t") if (d / "cells.tsv").exists() else None
        if cells is None:
            cnt = check_counts(cnt, smp, f"tet2/{d.name}")
        else:
            cnt = check_counts(cnt, cells, f"tet2/{d.name} cells", id_col="cell_id")
        omap = pd.read_csv(d / "orthologs.tsv", sep="\t") if (d / "orthologs.tsv").exists() else None
        out[d.name] = dict(samples=smp, counts=cnt, cells=cells, orthologs=omap)
    return out
