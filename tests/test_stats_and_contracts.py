import numpy as np
import pandas as pd
import pytest

from tet2isr import contracts, stats
from tet2isr.data_ingest.orthologs import map_counts, one_to_one
from tet2isr.isr_translation import dte
from tet2isr.tet2_transcriptome import analysis as tet2


def test_bh_monotone_and_bounded():
    p = np.array([0.001, 0.01, 0.02, 0.5, np.nan])
    q = stats.bh_fdr(p)
    assert np.isnan(q[4]) and (q[:4] >= p[:4]).all() and (np.diff(q[:4]) >= 0).all() and q[:4].max() <= 1


def test_lm_recovers_effect_and_rejects_unreplicated():
    rng = np.random.default_rng(0)
    g = np.r_[np.zeros(4), np.ones(4)]
    y = rng.normal(0, 0.1, (200, 8)) + 1.5 * g
    X = np.column_stack([np.ones(8), g])
    r = stats.lm_contrasts(y, X, {"g": np.array([0, 1.0])})["g"]
    assert abs(r.log2fc.mean() - 1.5) < 0.02 and ((r.ci_lo < 1.5) & (r.ci_hi > 1.5)).mean() > 0.9
    with pytest.raises(ValueError):
        stats.lm_contrasts(y[:, [0, 4]], X[[0, 4]], {"g": np.array([0, 1.0])})


def test_random_effects_meta_heterogeneity():
    est = np.array([[1.0, 1.0, 1.0], [0.0, 1.0, 2.0]])
    m = stats.random_effects_meta(est, np.full((2, 3), 0.1))
    assert m.tau[0] == 0 and m.tau[1] > 0.5 and m.pred_sd[1] > m.se[1] and m.I2[1] > 0.9


def test_ortholog_join_drops_ambiguous():
    o = pd.DataFrame(dict(source_gene_id=["a", "b", "b", "c", "d"], gene_id=["A", "B1", "B2", "C", "C"]))
    kept, info = one_to_one(o)
    assert list(kept.gene_id) == ["A"] and info["n_dropped"] == 4
    cnt = pd.DataFrame({"s1": [1, 2, 3, 4]}, index=["a", "b", "c", "d"])
    out, _ = map_counts(cnt, o)
    assert list(out.index) == ["A"] and out.s1.iloc[0] == 1


def test_count_sample_alignment_enforced(world):
    bad = world["cnt"].iloc[:, 1:]
    with pytest.raises(contracts.ContractError):
        contracts.check_counts(bad, world["smp"], "x")
    reordered = world["cnt"].iloc[:, ::-1]
    assert list(contracts.check_counts(reordered, world["smp"], "x").columns) == list(world["smp"].sample_id)


def test_rna_only_dataset_cannot_pass_as_paired(world):
    rna_only = world["smp"][world["smp"].assay == "rna"]
    with pytest.raises(contracts.ContractError):
        contracts.check_paired_rna_ribo(rna_only)


def test_dte_recovers_planted_truth(world):
    truth = pd.read_parquet(world["root"] / "truth/dte_by_context.parquet")
    for cid, d in world["sig"].groupby("context_id"):
        t = truth.loc[d.gene_id, cid].values
        assert np.corrcoef(d.dTE_log2fc, t)[0, 1] > 0.5
        # dRibo = dRNA + dTE exactly (contrast arithmetic)
        assert np.allclose(d.dRibo_log2fc, d.dRNA_log2fc + d.dTE_log2fc)


def test_pseudobulk_uses_samples_not_cells(world):
    s = contracts.load_tet2_studies(world["root"])["SYN_TET2_sc"]
    pb, meta = tet2.pseudobulk(s["counts"], s["cells"])
    assert pb.shape[1] == 12 and meta.sample_id.nunique() == 6          # 6 donors x 2 compartments
    assert pb.values.sum() == s["counts"].values.sum()
    de = tet2.genotype_de(pb[meta[meta.compartment == "HSC"].pb_id],
                          meta[meta.compartment == "HSC"].merge(s["samples"], on="sample_id").genotype)
    assert de.n_WT.iloc[0] == 3 and de.n_TET2.iloc[0] == 3 and de.df.iloc[0] < 50


def test_three_way_interaction_matches_difference_of_dtes(world):
    h = world["root"] / "tet2_isr_holdout"
    smp = pd.read_csv(h / "samples.tsv", sep="\t")
    o = dte.genotype_by_isr_dte(pd.read_parquet(h / "counts.parquet"), smp)
    assert np.allclose(o.I_TE, o.dTE_TET2 - o.dTE_WT)


def test_gene_space_mapping_sums_and_fills():
    from tet2isr.data_ingest.public import _to_gene_space
    df = pd.DataFrame({"s": [1, 2, 3, 4]}, index=["a", "b", "c", "zz"])
    out = _to_gene_space(df, pd.Series({"a": "G1", "b": "G1", "c": "G2"}), pd.Index(["G1", "G2", "G3"]))
    assert list(out.s) == [3, 3, 0]          # summed per gene, unmapped dropped, absent genes zero


def test_partial_spearman_removes_shared_technical_axis():
    from tet2isr.intersection import _partial_spearman
    rng = np.random.default_rng(1)
    z = pd.Series(rng.normal(size=3000))                      # shared covariate (e.g. gene length)
    x, y = z + rng.normal(size=3000), z + rng.normal(size=3000)
    raw, part, n = _partial_spearman(x, y, pd.DataFrame({"z": z}))
    assert raw > 0.4 and abs(part) < 0.06 and n == 3000
