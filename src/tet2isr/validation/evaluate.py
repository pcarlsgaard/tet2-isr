"""Held-out evaluation, baselines, ablations, leakage checks (brief section E)."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps

from ..modeling.models import VARIANTS, make_model

NOMINAL = (0.5, 0.8, 0.95)


def metrics(obs: np.ndarray, pred: np.ndarray, obs_fdr=None, sd=None, top_k: int = 50,
            programs: dict | None = None, index=None) -> dict:
    ok = np.isfinite(obs) & np.isfinite(pred)
    o, p = obs[ok] - np.median(obs[ok]), pred[ok] - np.median(pred[ok])
    const = np.allclose(p, p[0])
    out = dict(n_genes=int(ok.sum()),
               pearson=np.nan if const else float(np.corrcoef(o, p)[0, 1]),
               spearman=np.nan if const else float(sps.spearmanr(o, p)[0]),
               mae=float(np.abs(o - p).mean()), rmse=float(np.sqrt(((o - p) ** 2).mean())))
    sigm = (obs_fdr[ok] < 0.05) if obs_fdr is not None else np.abs(o) > 0.5
    out["n_sig_obs"] = int(sigm.sum())
    out["sign_concordance_sig"] = float((np.sign(o[sigm]) == np.sign(p[sigm])).mean()) if sigm.sum() and not const else np.nan
    k = min(top_k, len(o))
    to, tp = np.argsort(-np.abs(o))[:k], np.argsort(-np.abs(p))[:k]
    hit = np.intersect1d(to, tp)
    out[f"top{top_k}_precision"] = np.nan if const else float((np.sign(o[hit]) == np.sign(p[hit])).sum() / k)
    if programs and index is not None:
        idx = pd.Index(index[ok])
        po, pp = [], []
        for g in programs.values():
            m = idx.isin(g)
            if m.sum() >= 5:
                po.append(o[m].mean()); pp.append(p[m].mean())
        out["program_recovery_r"] = np.nan if const or len(po) < 3 else float(np.corrcoef(po, pp)[0, 1])
    if sd is not None:
        s = np.broadcast_to(sd, obs.shape)[ok]
        for lv in NOMINAL:
            zc = sps.norm.ppf(0.5 + lv / 2)
            out[f"coverage_{int(lv * 100)}"] = float((np.abs(o - p) <= zc * s).mean())
        out["mean_pi95_width"] = float((2 * 1.96 * s).mean())
    return out


def make_splits(ctx: pd.DataFrame, min_train_studies: int = 2) -> list[dict]:
    """Leakage rules: a held-out study removes ALL its contexts; a held-out stressor removes that
    stressor everywhere (dose variants share a stressor family); a held-out cell type removes
    every context of that cell type."""
    fam = ctx.stressor.str.replace(r"_(hi|lo|low|high)$", "", regex=True)
    out = []
    for kind, key in (("held_out_stressor", fam), ("held_out_cell_type", ctx.cell_type),
                      ("held_out_study", ctx.study_id)):
        for v in sorted(key.unique()):
            test = ctx.index[key == v]
            train = ctx.index[key != v]
            if ctx.loc[train, "study_id"].nunique() >= min_train_studies:
                out.append(dict(split=kind, held_out=v, train=list(train), test=list(test)))
    return out


def assert_no_leakage(split: dict, ctx: pd.DataFrame) -> None:
    tr, te = ctx.loc[split["train"]], ctx.loc[split["test"]]
    assert not set(tr.index) & set(te.index)
    if split["split"] == "held_out_study":
        assert not set(tr.study_id) & set(te.study_id), "study leakage"
    if split["split"] == "held_out_cell_type":
        assert not set(tr.cell_type) & set(te.cell_type), "cell-type leakage"
    if split["split"] == "held_out_stressor":
        f = lambda s: set(s.str.replace(r"_(hi|lo|low|high)$", "", regex=True))
        assert not f(tr.stressor) & f(te.stressor), "stressor leakage"


def run_heldout(sig, genes, ctx, state_cols, programs: pd.DataFrame, ref_genes, seed=0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    prog = {p: set(d.gene_id) for p, d in programs.groupby("program")}
    Y = sig.pivot(index="gene_id", columns="context_id", values="dTE_log2fc")
    F = sig.pivot(index="gene_id", columns="context_id", values="dTE_fdr")
    SE = sig.pivot(index="gene_id", columns="context_id", values="dTE_se")
    rows = []
    for sp in make_splits(ctx):
        assert_no_leakage(sp, ctx)
        tr_ctx, te_ctx = ctx.loc[sp["train"]], ctx.loc[sp["test"]]
        tr_sig = sig[sig.context_id.isin(sp["train"])]
        preds = {}
        for name in VARIANTS:
            m = make_model(name, genes, state_cols).fit(tr_sig, genes, tr_ctx)
            preds[name] = (m.predict(genes, te_ctx), m.fitted["resid_sd"])
        # null: gene labels permuted within every training context
        perm = tr_sig.copy()
        perm["dTE_log2fc"] = perm.groupby("context_id").dTE_log2fc.transform(lambda x: rng.permutation(x.values))
        m = make_model("M2_+cell_state", genes, state_cols).fit(perm, genes, tr_ctx)
        preds["NULL_permuted_labels"] = (m.predict(genes, te_ctx), m.fitted["resid_sd"])
        # baselines
        zero = pd.DataFrame(0.0, index=Y.index, columns=te_ctx.index)
        preds["B0_transcriptome_only/uniform_shift"] = (zero, None)
        Ytr = Y[sp["train"]]
        ztr = (tr_ctx[state_cols] - tr_ctx[state_cols].mean()) / tr_ctx[state_cols].std()
        zte = (te_ctx[state_cols] - tr_ctx[state_cols].mean()) / tr_ctx[state_cols].std()
        near = {c: Ytr[((ztr - zte.loc[c]) ** 2).sum(axis=1).idxmin()] for c in te_ctx.index}
        preds["B1_nearest_context"] = (pd.DataFrame(near), None)
        refm = Ytr.index.isin(ref_genes)
        eff = float((Ytr[refm].mean() - Ytr.median()).mean())
        preds["B2_reference_gene_set"] = (pd.DataFrame({c: np.where(refm, eff, 0.0) for c in te_ctx.index},
                                                       index=Ytr.index), None)
        for c in te_ctx.index:
            o = Y[c].dropna()
            for name, (P, rsd) in preds.items():
                p = P[c].reindex(o.index).values
                sd = None if rsd is None else np.sqrt(rsd ** 2 + 0 * SE[c].reindex(o.index).values)
                rows.append(dict(split=sp["split"], held_out=sp["held_out"], test_context=c,
                                 stressor_class=te_ctx.loc[c, "stressor_class"],
                                 cell_type_in_training=bool(te_ctx.loc[c, "cell_type"] in set(tr_ctx.cell_type)),
                                 n_train_contexts=len(tr_ctx), model=name,
                                 **metrics(o.values, p, F[c].reindex(o.index).values, sd,
                                           programs=prog, index=o.index.values)))
    return pd.DataFrame(rows)


def gene_family_leakage(sig, genes, ctx, n_folds=5, seed=0) -> pd.DataFrame:
    """Sequence-only prediction of the pooled unit response for UNSEEN genes, with folds drawn
    (a) at random over genes [paralogs can straddle train/test] and (b) grouped by gene family."""
    rng = np.random.default_rng(seed)
    d = sig.drop(columns=["intensity"], errors="ignore").merge(ctx[["intensity"]], left_on="context_id", right_index=True)
    y = (d.assign(u=d.dTE_log2fc / d.intensity).groupby("gene_id").u.mean())
    g = genes.set_index("gene_id").loc[y.index]
    feats = [c for c in g.columns if c.startswith("seq_")]
    X = ((g[feats] - g[feats].mean()) / g[feats].std()).values
    X = np.column_stack([np.ones(len(X)), X])
    rows = []
    for scheme in ("random_gene_split", "family_grouped_split"):
        keys = np.arange(len(g)) if scheme == "random_gene_split" else pd.factorize(g.family_id)[0]
        fold = rng.integers(0, n_folds, keys.max() + 1)[keys]
        pred = np.empty(len(y))
        for f in range(n_folds):
            tr, te = fold != f, fold == f
            b = np.linalg.solve(X[tr].T @ X[tr] + np.eye(X.shape[1]), X[tr].T @ y.values[tr])
            pred[te] = X[te] @ b
        fam_in_both = 0 if scheme != "random_gene_split" else int(
            sum(len(set(fold[g.family_id.values == fam])) > 1 for fam in g.family_id[g.family_id.duplicated()].unique()))
        rows.append(dict(scheme=scheme, pearson=float(np.corrcoef(pred, y.values)[0, 1]),
                         rmse=float(np.sqrt(((pred - y.values) ** 2).mean())),
                         families_split_across_folds=fam_in_both))
    # nearest-paralog memorisation: how much could a model gain purely from family identity?
    fam_mean = y.groupby(g.family_id.values).transform(lambda x: (x.sum() - x) / max(len(x) - 1, 1) if len(x) > 1 else np.nan)
    ok = fam_mean.notna()
    rows.append(dict(scheme="paralog_copy_upper_bound(genes with paralogs only)",
                     pearson=float(np.corrcoef(fam_mean[ok], y[ok])[0, 1]),
                     rmse=float(np.sqrt(((fam_mean[ok] - y[ok]) ** 2).mean())),
                     families_split_across_folds=int(ok.sum())))
    return pd.DataFrame(rows)


def factor_stability(sig, genes, ctx, state_cols, rank=2, n_boot=30, seed=0) -> pd.DataFrame:
    """|cosine| between low-rank gene loadings from the full fit and context-bootstrap refits."""
    rng = np.random.default_rng(seed)
    ref = make_model("M3_low_rank", genes, state_cols).fit(sig, genes, ctx).fitted
    if "L" not in ref or len(ctx) < rank + 4:
        return pd.DataFrame([dict(factor="not fitted", note=f"{len(ctx)} contexts: too few for a rank-{rank} residual model")])
    L0 = ref["L"] / np.linalg.norm(ref["L"].values, axis=0)
    rows = []
    for b in range(n_boot):
        keep = rng.choice(ctx.index, len(ctx) - 2, replace=False)  # delete-2 jackknife over contexts
        f = make_model("M3_low_rank", genes, state_cols).fit(sig[sig.context_id.isin(keep)], genes, ctx.loc[keep]).fitted
        L = (f["L"] / np.linalg.norm(f["L"].values, axis=0)).reindex(L0.index).fillna(0)
        cs = np.abs(L0.values.T @ L.values)
        for r in range(rank):
            rows.append(dict(replicate=b, factor=r + 1, max_abs_cosine=float(cs[r].max())))
    sv = ref["sv"]
    out = pd.DataFrame(rows).groupby("factor").max_abs_cosine.agg(["mean", "min"]).reset_index()
    out["variance_fraction_of_residual"] = (sv[:rank] ** 2 / (sv ** 2).sum())
    return out


def genotype_holdout(obs: pd.DataFrame, m1_gene: pd.DataFrame, m2_gene: pd.DataFrame) -> pd.DataFrame:
    """Compare predictions with a directly measured genotype x ISR translatome (never used in fitting)."""
    o = obs.set_index("gene_id")
    rows = []
    for (model, d) in (("M1", m1_gene), ("M2", m2_gene)):
        for (sid, inten), dd in d.groupby(["tet2_study", "isr_intensity"]):
            dd = dd.set_index("gene_id").reindex(o.index)
            if model == "M1":
                tw = tk = dd.dTE_transferred.values
                ite = np.zeros(len(o))
                sd = None
            else:
                tw, tk, ite = dd.dTE_WT.values, dd.dTE_TET2.values, dd.I_TE.values
                sd = (dd.I_TE_hi - dd.I_TE_lo).values / (2 * 1.96)
            for target, obsv, fdr, pred, s in (("dTE_WT", o.dTE_WT, o.dTE_WT_fdr, tw, None),
                                               ("dTE_TET2", o.dTE_TET2, o.dTE_TET2_fdr, tk, None),
                                               ("I_TE", o.I_TE, o.I_TE_fdr, ite,
                                                None if sd is None else np.sqrt(sd ** 2 + o.I_TE_se.values ** 2))):
                rows.append(dict(model=model, tet2_study=sid, isr_intensity=inten, target=target,
                                 **metrics(obsv.values, pred, fdr.values, s)))
    return pd.DataFrame(rows)


def holdout_program_level(obs: pd.DataFrame, m2_gene: pd.DataFrame, programs: pd.DataFrame,
                          intensity_label: str = "moderate") -> pd.DataFrame:
    """Gene-set-level comparison of observed vs predicted I_TE. Per-gene I_TE is usually
    noise-limited; set means are the realistic unit of validation.  SE treats genes as
    exchangeable (ignores gene-gene correlation) and is therefore optimistic."""
    o = obs.set_index("gene_id")
    rows = []
    for p, d in programs.groupby("program"):
        g = o.index.intersection(d.gene_id)
        if len(g) < 5:
            continue
        x = o.loc[g, "I_TE"] - o.I_TE.median()
        base = dict(program=p, n_genes=len(g), observed_I_TE=float(x.mean()),
                    observed_se=float(x.std(ddof=1) / np.sqrt(len(g))))
        for (src, sid), dd in m2_gene[m2_gene.isr_intensity == intensity_label].groupby(
                ["intensity_source", "tet2_study"]):
            pr = dd.set_index("gene_id").reindex(g)
            rows.append(dict(**base, intensity_source=src, tet2_study=sid,
                             predicted_I_TE_M2=float(pr.I_TE.mean()),
                             predicted_lo=float(pr.I_TE_lo.mean()), predicted_hi=float(pr.I_TE_hi.mean()),
                             predicted_I_TE_M1=0.0))
    return pd.DataFrame(rows)


def unengaged_check(sig_all, sig, genes, ctx, state_cols, programs) -> pd.DataFrame:
    """Negative control: how well does the signature learned from ISR-engaged contexts 'predict' dTE
    in stress contexts that failed the engagement gate?  High agreement would mean the signature is
    generic stress/batch structure rather than ISR-specific."""
    rest = sig_all[~sig_all.context_id.isin(ctx.index)]
    if rest.empty:
        return pd.DataFrame()
    m = make_model("M1b_signature+dose", genes, state_cols).fit(sig, genes, ctx)
    u = m.fitted["u"]
    prog = {p: set(d.gene_id) for p, d in programs.groupby("program")}
    rows = []
    for cid, d in rest.groupby("context_id"):
        o = d.set_index("gene_id").dTE_log2fc
        g = o.index.intersection(u.index)
        rows.append(dict(context_id=cid, **metrics(o[g].values, u[g].values, d.set_index("gene_id").dTE_fdr[g].values,
                                                   programs=prog, index=g.values)))
    return pd.DataFrame(rows)
