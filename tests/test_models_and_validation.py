import numpy as np
import pandas as pd

from tet2isr import contracts
from tet2isr.isr_translation import dte
from tet2isr.modeling import models
from tet2isr.tet2_transcriptome import analysis as tet2
from tet2isr.validation import evaluate

from conftest import STATE


def test_splits_have_no_leakage(world):
    splits = evaluate.make_splits(world["ctx"])
    assert {s["split"] for s in splits} == {"held_out_stressor", "held_out_cell_type", "held_out_study"}
    for s in splits:
        evaluate.assert_no_leakage(s, world["ctx"])
    # dose variants of one stressor are held out together
    tg = [s for s in splits if s["split"] == "held_out_stressor" and s["held_out"] == "thapsigargin"][0]
    assert not any("thapsigargin" in c for c in tg["train"])


def test_model_shapes_and_unseen_gene_fallback(world):
    g, c, sig = world["genes"], world["ctx"], world["sig"]
    for name in models.VARIANTS:
        m = models.make_model(name, g, STATE).fit(sig, g, c)
        P = m.predict(g, c)
        assert P.shape == (len(g), len(c)) and np.isfinite(P.values).all()
    m = models.make_model("M3_low_rank", g, STATE).fit(sig, g, c)
    assert m.fitted["L"].shape[1] == 2
    seen = sig[sig.gene_id.isin(g.gene_id[:300])]
    m = models.make_model("M1b_signature+dose", g, STATE).fit(seen, g, c)
    assert np.isfinite(m.predict(g, c).values).all()        # unseen genes fall back to sequence prior


def test_m1_has_zero_te_interaction_by_construction(world):
    r = tet2.run_all(contracts.load_tet2_studies(world["root"]))
    pooled = dte.pooled_signature(world["sig"], per_unit_intensity=True)
    gene, prog = models.m1_predict(r["expr"], pooled, world["programs"], {"moderate": 1.0}, n_mc=20)
    assert (gene.I_TE == 0).all() and (gene.I_S_log2_gene == 0).all()
    assert np.allclose(gene.I_S_linear, gene.dS_linear_TET2 - gene.dS_linear_WT)
    assert (gene.evidence == "transferred").all()
    assert (prog.I_S_log2_program_lo <= prog.I_S_log2_program_hi).all()


def test_m2_predicts_no_interaction_without_state_difference(world):
    """If WT and TET2 have the same basal state, M2 must return I_TE == 0."""
    g, c, sig = world["genes"], world["ctx"], world["sig"]
    m = models.make_model("M2_+cell_state", g, STATE).fit(sig, g, c)
    z = pd.DataFrame([c[STATE].mean(), c[STATE].mean()], index=["WT", "TET2"])
    z["stressor_class"] = "canonical"
    u = m.unit_response(g, z)
    assert np.allclose(u["WT"], u["TET2"])
    far = z.copy()
    far.loc["TET2", STATE[0]] = c[STATE[0]].max() + 5
    assert list(m.extrapolation(far)) == [False, True]


def test_permuted_labels_give_no_signal(world):
    ref = world["programs"].loc[world["programs"].program == "isr_uorf_targets", "gene_id"]
    ev = evaluate.run_heldout(world["sig"], world["genes"], world["ctx"], STATE, world["programs"], ref)
    s = ev[ev.split == "held_out_study"].groupby("model").pearson.mean()
    assert abs(s["NULL_permuted_labels"]) < 0.15
    assert s["M1b_signature+dose"] > 0.4 > s["NULL_permuted_labels"]


def test_family_grouped_split_separates_paralogs(world):
    r = evaluate.gene_family_leakage(world["sig"], world["genes"], world["ctx"]).set_index("scheme")
    assert r.loc["family_grouped_split", "families_split_across_folds"] == 0
    assert r.loc["random_gene_split", "families_split_across_folds"] > 0
