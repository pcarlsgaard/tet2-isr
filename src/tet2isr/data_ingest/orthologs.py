"""One-to-one ortholog mapping. Ambiguous (one-to-many / many-to-one) pairs are dropped,
never summed, and the number dropped is reported."""
from __future__ import annotations

import pandas as pd


def one_to_one(omap: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    omap = omap[["source_gene_id", "gene_id"]].dropna().drop_duplicates()
    ok = (~omap.source_gene_id.duplicated(keep=False)) & (~omap.gene_id.duplicated(keep=False))
    kept = omap[ok]
    return kept, dict(n_pairs=len(omap), n_one_to_one=len(kept), n_dropped=int((~ok).sum()))


def map_counts(counts: pd.DataFrame, omap: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    kept, info = one_to_one(omap)
    kept = kept[kept.source_gene_id.isin(counts.index)]
    out = counts.loc[kept.source_gene_id.values]
    out.index = pd.Index(kept.gene_id.values, name="gene_id")
    info["n_mapped_genes"] = len(out)
    return out, info
