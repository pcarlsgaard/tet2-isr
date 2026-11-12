import numpy as np

from tet2isr import effects as fx


def test_shared_dte_gives_zero_interactions_but_nonzero_linear():
    m_wt, m_ko, te0, te1 = 10.0, 20.0, 1.0, 2.0
    s = lambda m, te: fx.synthesis(m, te)
    assert fx.interaction_te(te0, te1, te0, te1) == 0
    assert np.isclose(fx.interaction_log_s(s(m_wt, te0), s(m_wt, te1), s(m_ko, te0), s(m_ko, te1)), 0)
    assert fx.interaction_linear_s(s(m_wt, te0), s(m_wt, te1), s(m_ko, te0), s(m_ko, te1)) == 10.0


def test_true_te_interaction_is_detected():
    assert np.isclose(fx.interaction_te(1, 2, 1, 4), 1.0)
    assert np.isclose(fx.interaction_log_s(10, 20, 20, 80), 1.0)


def test_aggregate_interaction_arises_from_composition_alone():
    d = np.array([1.0, -1.0])
    wt, ko = np.array([10.0, 40.0]), np.array([20.0, 40.0])
    assert not np.isclose(fx.aggregate_log_change(wt, d), fx.aggregate_log_change(ko, d))
    # identical composition (scaled) -> identical aggregate change
    assert np.isclose(fx.aggregate_log_change(wt, d), fx.aggregate_log_change(3 * wt, d))


def test_two_by_two_table():
    s = fx.two_by_two_example()["summary"].set_index("unit")
    assert s.loc["A", "I_TE"] == 0 and s.loc["A", "I_S_log2"] == 0 and s.loc["A", "I_S_linear"] == 10
    assert s.loc["B", "I_S_linear"] == 0
    assert s.loc["pathway{A,B}", "I_S_log2"] > 0.3
