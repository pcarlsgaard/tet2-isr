"""Model tiers M0-M3 (brief section D).

M0  transcriptome only (no translation layer).
M1  fixed transferred dTE applied to genotype-specific mRNA pools. I_TE == 0 BY ASSUMPTION.
M2  dTE[g,c] = I_c * ( u_g + sum_j z_cj * [1, X_g] @ gamma_j ),  ridge-regularised;
    u_g = gene-level unit response (or X_g @ beta for unseen genes), z_c = basal cell-state
    program scores of context c, X_g = sequence (optionally epigenomic) features.
M3  low-rank: residual gene x context matrix factorised by SVD; context scores regressed on z_c.
M4  intentionally not implemented (see docs/decisions.md).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from ..effects import aggregate_log_change
from ..stats import log_cpm, size_factors


# ------------------------------------------------------------------ shared helpers
def program_scores(log2cpm: pd.DataFrame, programs: pd.DataFrame, names) -> pd.DataFrame:
    """Per-sample program score: mean log2CPM of members minus the sample's all-gene mean."""
    out = {}
    for p in names:
        g = programs.loc[programs.program == p, "gene_id"]
        g = g[g.isin(log2cpm.index)]
        out[p] = log2cpm.loc[g].mean(axis=0) - log2cpm.mean(axis=0)
    return pd.DataFrame(out)


def context_states(counts, samples, programs, names) -> pd.DataFrame:
    """Basal state z_c of every ISR context from its CONTROL RNA libraries only."""
    s = samples[(samples.assay == "rna") & (samples.condition == "control")]
    lc = pd.DataFrame(log_cpm(counts[s.sample_id].values), index=counts.index, columns=s.sample_id)
    sc = program_scores(lc, programs, names)
    sc["context_id"] = s.context_id.values
    return sc.groupby("context_id").mean()


def estimate_intensity(sig: pd.DataFrame, ref_genes) -> pd.Series:
    """Data-driven ISR intensity per context: mean dTE of reference uORF targets minus the
    context median, scaled to mean 1.  Doses of different stressors are not commensurable."""
    r = sig[sig.gene_id.isin(ref_genes)].groupby("context_id").dTE_log2fc.mean()
    r = r - sig.groupby("context_id").dTE_log2fc.median()
    return (r / r.mean()).clip(lower=0.05).rename("intensity")


def isr_engagement(sig: pd.DataFrame, ref_genes, min_dte: float = 0.5, alpha: float = 0.05) -> pd.DataFrame:
    """Positive-control gate. A context counts as ISR-engaged only if the reference uORF-target set
    has mean dTE (relative to the context median) >= min_dte with one-sided t-test p < alpha across
    the reference genes. Unengaged contexts are reported but do not define the ISR signature."""
    from scipy import stats as sps
    med = sig.groupby("context_id").dTE_log2fc.median()
    rows = []
    for cid, d in sig[sig.gene_id.isin(ref_genes)].groupby("context_id"):
        x = d.dTE_log2fc.values - med[cid]
        p = sps.ttest_1samp(x, 0, alternative="greater").pvalue if len(x) > 2 else np.nan
        rows.append(dict(context_id=cid, n_ref_genes=len(x), ref_mean_dTE=float(x.mean()),
                         ref_frac_positive=float((x > 0).mean()), ref_p_one_sided=float(p),
                         engaged=bool(x.mean() >= min_dte and p < alpha)))
    return pd.DataFrame(rows).set_index("context_id")


# ------------------------------------------------------------------ M0
def m0_programs(expr: dict, programs: pd.DataFrame) -> pd.DataFrame:
    """Replicate-aware program-level genotype contrasts (Welch t on per-sample scores)."""
    from scipy import stats as sps
    rows = []
    for sid, e in expr.items():
        c = e["counts"]
        lc = pd.DataFrame(log_cpm(c.values, size_factors(c.values)), index=c.index, columns=c.columns)
        sc = program_scores(lc, programs, programs.program.unique())
        t = (e["genotype"] == "TET2").values
        for p in sc.columns:
            a, b = sc[p].values[t], sc[p].values[~t]
            d = a.mean() - b.mean()
            se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
            rows.append(dict(study_id=sid, program=p, log2fc=d, se=se,
                             ci_lo=d - 1.96 * se, ci_hi=d + 1.96 * se,
                             pvalue=sps.ttest_ind(a, b, equal_var=False).pvalue,
                             n_WT=int((~t).sum()), n_TET2=int(t.sum()), evidence="observed"))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ M1
def m1_predict(expr: dict, pooled: pd.DataFrame, programs: pd.DataFrame, intensities: dict,
               n_mc: int = 300, seed: int = 0):
    """Overlay a common dTE (per unit intensity; pooled across training contexts, SD = pred_sd
    including between-context heterogeneity) on WT and TET2 basal mRNA pools.

    Early-phase assumption: mRNA abundance is unchanged by stress; basal TE is set to 1 for
    every gene (arbitrary common scale), so linear outputs are 'CPM x relative TE' units, not
    molecules.  Gene-level I_TE and gene-level log-scale I_S are zero by construction.
    """
    rng = np.random.default_rng(seed)
    pooled = pooled.set_index("gene_id")
    gene_rows, prog_rows = [], []
    for sid, e in expr.items():
        c = e["counts"]
        cpm = pd.DataFrame(2 ** log_cpm(c.values, size_factors(c.values)), index=c.index, columns=c.columns)
        genes = cpm.index.intersection(pooled.index)
        cpm = cpm.loc[genes]
        wt, ko = cpm.loc[:, (e["genotype"] == "WT").values], cpm.loc[:, (e["genotype"] == "TET2").values]
        mu, sd = pooled.loc[genes, "dTE"].values, pooled.loc[genes, "pred_sd"].values
        for label, inten in intensities.items():
            d = inten * mu
            mw, mk = wt.mean(axis=1).values, ko.mean(axis=1).values
            draws = np.empty((n_mc, len(genes)))
            pdraw = {p: np.empty((n_mc, 4)) for p in programs.program.unique()}
            pidx = {p: np.isin(genes, programs.loc[programs.program == p, "gene_id"]) for p in pdraw}
            for b in range(n_mc):
                db = inten * rng.normal(mu, sd)
                mwb = wt.iloc[:, rng.integers(0, wt.shape[1], wt.shape[1])].mean(axis=1).values
                mkb = ko.iloc[:, rng.integers(0, ko.shape[1], ko.shape[1])].mean(axis=1).values
                draws[b] = (mkb - mwb) * (2 ** db - 1)
                tot_w, tot_k = (mwb * 2 ** db).sum() / mwb.sum(), (mkb * 2 ** db).sum() / mkb.sum()
                for p, ix in pidx.items():
                    lw, lk = aggregate_log_change(mwb[ix], db[ix]), aggregate_log_change(mkb[ix], db[ix])
                    # share of total synthesis potential, log2 change
                    pdraw[p][b] = (lw, lk, lk - lw, (lk - np.log2(tot_k)) - (lw - np.log2(tot_w)))
            lo, hi = np.percentile(draws, [2.5, 97.5], axis=0)
            gene_rows.append(pd.DataFrame(dict(
                gene_id=genes, tet2_study=sid, isr_intensity=label, model="M1",
                M_WT_cpm=mw, M_TET2_cpm=mk, basal_log2fc_tet2_vs_wt=np.log2(mk / mw),
                dTE_transferred=d, dTE_pred_sd=inten * sd,
                dS_linear_WT=mw * (2 ** d - 1), dS_linear_TET2=mk * (2 ** d - 1),
                I_S_linear=(mk - mw) * (2 ** d - 1), I_S_linear_lo=lo, I_S_linear_hi=hi,
                I_S_log2_gene=0.0, I_TE=0.0,
                I_TE_note="zero by M1 assumption, not an estimate",
                evidence="transferred")))
            for p, ix in pidx.items():
                q = np.percentile(pdraw[p], [2.5, 50, 97.5], axis=0)
                prog_rows.append(dict(
                    program=p, tet2_study=sid, isr_intensity=label, model="M1", n_genes=int(ix.sum()),
                    dlog2S_WT=q[1, 0], dlog2S_TET2=q[1, 1],
                    I_S_log2_program=q[1, 2], I_S_log2_program_lo=q[0, 2], I_S_log2_program_hi=q[2, 2],
                    I_S_share_log2=q[1, 3], I_S_share_lo=q[0, 3], I_S_share_hi=q[2, 3],
                    mechanism="substrate/compositional only (abundance-weighted aggregation of a common dTE)",
                    evidence="transferred"))
    return pd.concat(gene_rows, ignore_index=True), pd.DataFrame(prog_rows)


# ------------------------------------------------------------------ M2 / M3
@dataclass
class TransferModel:
    """Ridge transfer model. Flags select ablations."""
    feature_cols: list
    state_cols: list
    use_gene_effect: bool = True       # False -> sequence-only (works for unseen genes)
    use_intensity: bool = True         # False -> every context treated as intensity 1
    use_state: bool = True             # basal cell-state x sequence interactions
    use_stress_class: bool = False     # pleiotropic-vs-canonical x sequence interaction
    low_rank: int = 0                  # >0 -> M3 residual SVD of this rank
    alpha_seq: float = 1.0
    alpha_int: float = 30.0
    shrink_gene: float = 0.5           # pseudo-observations shrinking u_g to the sequence prediction
    fitted: dict = field(default_factory=dict)

    def _X(self, genes: pd.DataFrame):
        X = genes.set_index("gene_id")[self.feature_cols].astype(float)
        mu, sd = self.fitted.get("x_mu"), self.fitted.get("x_sd")
        if mu is None:
            mu, sd = X.mean(), X.std().replace(0, 1)
            self.fitted.update(x_mu=mu, x_sd=sd)
        X = (X - mu) / sd
        X.insert(0, "const", 1.0)
        return X

    def _z(self, ctx: pd.DataFrame):
        cols = list(self.state_cols) if self.use_state else []
        Z = ctx[cols].astype(float).copy()
        if "z_mu" not in self.fitted:
            self.fitted.update(z_mu=Z.mean(), z_sd=Z.std().replace(0, 1))
        Z = (Z - self.fitted["z_mu"]) / self.fitted["z_sd"]
        if self.use_stress_class:
            Z["pleiotropic"] = (ctx.stressor_class == "pleiotropic").astype(float).values
        return Z

    def fit(self, sig: pd.DataFrame, genes: pd.DataFrame, ctx: pd.DataFrame):
        """sig: gene_id, context_id, dTE_log2fc.  ctx: indexed by context_id with intensity,
        state columns and stressor_class."""
        self.fitted = {}
        X = self._X(genes)
        Z = self._z(ctx)
        inten = ctx.intensity if self.use_intensity else pd.Series(1.0, index=ctx.index)
        Y = sig.pivot(index="gene_id", columns="context_id", values="dTE_log2fc")
        Y = Y.reindex(index=X.index.intersection(Y.index), columns=ctx.index)
        Y = Y - Y.median(axis=0)                       # TE is relative: centre each context
        Xg = X.loc[Y.index].values
        I = inten.loc[Y.columns].values
        Zc = Z.loc[Y.columns].values
        mask = np.isfinite(Y.values)
        Y0 = np.where(mask, Y.values, 0.0)
        p, q = Xg.shape[1], Zc.shape[1]
        gamma = np.zeros((q, p))
        u = np.zeros(len(Y))
        beta = np.zeros(p)
        ridge = lambda A, b, a: np.linalg.solve(A.T @ A + a * np.eye(A.shape[1]), A.T @ b)
        for _ in range(4 if q else 1):
            R = Y0 - I[None, :] * (Xg @ gamma.T @ Zc.T) if q else Y0
            num = (mask * R * I[None, :]).sum(axis=1)
            den = (mask * I[None, :] ** 2).sum(axis=1)
            beta = ridge(Xg, num / np.maximum(den, 1e-9), self.alpha_seq)
            prior = Xg @ beta
            k = self.shrink_gene
            u = (num + k * prior) / (den + k) if self.use_gene_effect else prior
            if q:
                R2 = (Y0 - I[None, :] * u[:, None]) * mask
                # design for gamma: rows (g,c), cols (j,f) = I_c z_cj X_gf  -> solve via normal eqs
                W = (I[:, None] * Zc)                                   # C x q
                A = np.einsum("gc,cj,ck,gf,gh->jfkh", mask.astype(float), W, W, Xg, Xg).reshape(q * p, q * p)
                b = np.einsum("gc,cj,gf->jf", R2, W, Xg).reshape(q * p)
                gamma = np.linalg.solve(A + self.alpha_int * len(Y) / 100 * np.eye(q * p), b).reshape(q, p)
        resid = (Y0 - I[None, :] * (u[:, None] + (Xg @ gamma.T @ Zc.T if q else 0))) * mask
        L = None
        if self.low_rank and len(Y.columns) > self.low_rank + 1:
            Rn = resid / I[None, :]
            U, S, Vt = np.linalg.svd(Rn, full_matrices=False)
            r = self.low_rank
            L, sc = U[:, :r] * S[:r], Vt[:r].T                          # G x r, C x r
            Zs = self._z_all(ctx).loc[Y.columns].values
            Zd = np.column_stack([np.ones(len(Zs)), Zs])
            B = np.linalg.solve(Zd.T @ Zd + 1.0 * np.eye(Zd.shape[1]), Zd.T @ sc)
            self.fitted.update(L=pd.DataFrame(L, index=Y.index), B=B, sv=S)
            resid = resid - I[None, :] * (L @ (Zd @ B).T) * mask
        self.fitted.update(u=pd.Series(u, index=Y.index), beta=beta, gamma=gamma,
                           resid_sd=float(np.sqrt((resid[mask] ** 2).mean())),
                           mean_intensity=float(I.mean()),
                           z_range=(Z.min(), Z.max()), train_contexts=list(Y.columns))
        return self

    def _z_all(self, ctx):
        Z = ctx[list(self.state_cols)].astype(float)
        mu = self.fitted.setdefault("zall_mu", Z.mean())
        sd = self.fitted.setdefault("zall_sd", Z.std().replace(0, 1))
        return (Z - mu) / sd

    def unit_response(self, genes: pd.DataFrame, ctx: pd.DataFrame) -> pd.DataFrame:
        """Predicted dTE per unit intensity, genes x contexts."""
        X = self._X(genes)
        Z = self._z(ctx)
        u = self.fitted["u"].reindex(X.index)
        prior = pd.Series(X.values @ self.fitted["beta"], index=X.index)
        u = u.fillna(prior)                                              # unseen genes
        out = np.repeat(u.values[:, None], len(ctx), axis=1)
        if Z.shape[1]:
            out = out + X.values @ self.fitted["gamma"].T @ Z.values.T
        if "L" in self.fitted:
            Zd = np.column_stack([np.ones(len(ctx)), self._z_all(ctx).values])
            L = self.fitted["L"].reindex(X.index).fillna(0.0).values
            out = out + L @ (Zd @ self.fitted["B"]).T
        return pd.DataFrame(out, index=X.index, columns=ctx.index)

    def predict(self, genes, ctx, intensity=None) -> pd.DataFrame:
        if intensity is None:
            intensity = np.full(len(ctx), self.fitted["mean_intensity"] if self.use_intensity else 1.0)
        return self.unit_response(genes, ctx) * np.asarray(intensity)[None, :]

    def extrapolation(self, ctx: pd.DataFrame) -> pd.Series:
        """True where any state covariate lies outside the range of the training contexts."""
        Z = self._z(ctx)
        lo, hi = self.fitted["z_range"]
        cols = [c for c in Z.columns if c in lo.index]
        return ((Z[cols] < lo[cols]) | (Z[cols] > hi[cols])).any(axis=1)


VARIANTS = {  # name -> kwargs; ordered as the ablation ladder
    "M2a_sequence_only": dict(use_gene_effect=False, use_intensity=True, use_state=False),
    "M1_gene_signature": dict(use_gene_effect=True, use_intensity=False, use_state=False),
    "M1b_signature+dose": dict(use_gene_effect=True, use_intensity=True, use_state=False),
    "M2_+cell_state": dict(use_gene_effect=True, use_intensity=True, use_state=True),
    "M2c_+stress_class": dict(use_gene_effect=True, use_intensity=True, use_state=True, use_stress_class=True),
    "M2d_+epigenome": dict(use_gene_effect=True, use_intensity=True, use_state=True, _epi=True),
    "M3_low_rank": dict(use_gene_effect=True, use_intensity=True, use_state=False, low_rank=2),
}


def make_model(name: str, genes: pd.DataFrame, state_cols) -> TransferModel:
    kw = dict(VARIANTS[name])
    epi = kw.pop("_epi", False)
    feats = [c for c in genes.columns if c.startswith("seq_") or (epi and c.startswith("epi_"))]
    return TransferModel(feature_cols=feats, state_cols=list(state_cols), **kw)


def bootstrap_fit(name, sig, genes, ctx, state_cols, n_boot=30, seed=0):
    """Refit on resampled SOURCE STUDIES (not genes, not libraries)."""
    rng = np.random.default_rng(seed)
    studies = ctx.study_id.unique()
    fits = []
    for _ in range(n_boot):
        pick = rng.choice(studies, len(studies), replace=True)
        if len(set(pick)) < min(3, len(studies)):
            continue
        parts, cparts = [], []
        for i, s in enumerate(pick):  # duplicated studies get distinct context ids
            c = ctx[ctx.study_id == s]
            ren = {k: f"{k}#b{i}" for k in c.index}
            cparts.append(c.rename(index=ren))
            parts.append(sig[sig.context_id.isin(c.index)].assign(
                context_id=lambda d: d.context_id.map(ren)))
        fits.append(make_model(name, genes, state_cols).fit(pd.concat(parts), genes, pd.concat(cparts)))
    return fits


def m2_genotype_prediction(sig, genes, ctx, state_cols, expr: dict, programs, intensities: dict,
                           n_boot=40, seed=0, name="M2_+cell_state"):
    """Model-based I_TE: the part of the learned state-dependence of dTE that is engaged by the
    WT-vs-TET2 difference in basal state.  Unvalidated unless a direct translatome exists."""
    rng = np.random.default_rng(seed)
    fits = bootstrap_fit(name, sig, genes, ctx, state_cols, n_boot, seed)
    full = make_model(name, genes, state_cols).fit(sig, genes, ctx)
    rows, states = [], []
    for sid, e in expr.items():
        c = e["counts"]
        lc = pd.DataFrame(log_cpm(c.values, size_factors(c.values)), index=c.index, columns=c.columns)
        sc = program_scores(lc, programs, state_cols)
        t = (e["genotype"] == "TET2").values

        def geno_ctx(scw, sck):
            z = pd.DataFrame([scw.mean(), sck.mean()], index=["WT", "TET2"])
            z["stressor_class"] = "canonical"
            return z
        z_full = geno_ctx(sc[~t], sc[t])
        flag = full.extrapolation(z_full)
        states.append(z_full[state_cols].assign(tet2_study=sid, outside_training_range=flag.values))
        ur = full.unit_response(genes, z_full)
        draws = []
        for f in fits:  # joint bootstrap: ISR source studies x TET2 biological replicates
            w, k = sc[~t], sc[t]
            zb = geno_ctx(w.iloc[rng.integers(0, len(w), len(w))], k.iloc[rng.integers(0, len(k), len(k))])
            u = f.unit_response(genes, zb)
            draws.append(np.column_stack([u["WT"], u["TET2"], u["TET2"] - u["WT"]]))
        D = np.stack(draws)                                        # B x G x 3
        ite = D[:, :, 2]
        rs = pd.DataFrame(ite.T).corr(method="spearman").values
        rank_stab = float(rs[np.triu_indices_from(rs, 1)].mean())
        point = (ur["TET2"] - ur["WT"]).values
        sign_stab = (np.sign(ite) == np.sign(point)[None, :]).mean(axis=0)
        for label, inten in intensities.items():
            q = np.percentile(D * inten, [2.5, 97.5], axis=0)
            rows.append(pd.DataFrame(dict(
                gene_id=ur.index, tet2_study=sid, isr_intensity=label, model=name,
                dTE_WT=inten * ur["WT"].values, dTE_WT_lo=q[0, :, 0], dTE_WT_hi=q[1, :, 0],
                dTE_TET2=inten * ur["TET2"].values, dTE_TET2_lo=q[0, :, 1], dTE_TET2_hi=q[1, :, 1],
                I_TE=inten * point, I_TE_lo=q[0, :, 2], I_TE_hi=q[1, :, 2],
                sign_stability=sign_stab, rank_stability_spearman=rank_stab,
                extrapolation_flag=bool(flag.any()),
                evidence="transferred", I_TE_note="model-based; not validated by direct measurement")))
    return pd.concat(rows, ignore_index=True), pd.concat(states), full
