"""Small, dependency-light statistics used across the pipeline.

Count models here are log-CPM linear models with limma-style empirical-Bayes
variance moderation.  They are adequate for the smoke test and for a first pass
on processed count tables; DESeq2 ``deltaTE`` / anota2seq remain the recommended
confirmatory analysis (see docs/methods.md).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import special, stats as sps


def size_factors(counts: np.ndarray) -> np.ndarray:
    """Median-of-ratios size factors. counts: genes x samples."""
    c = np.asarray(counts, float)
    ok = (c > 0).all(axis=1)
    if ok.sum() < 10:
        lib = c.sum(axis=0)
        return lib / np.exp(np.mean(np.log(lib)))
    lg = np.log(c[ok])
    return np.exp(np.median(lg - lg.mean(axis=1, keepdims=True), axis=0))


def log_cpm(counts: np.ndarray, sf: np.ndarray | None = None, prior: float = 0.5) -> np.ndarray:
    c = np.asarray(counts, float)
    sf = size_factors(c) if sf is None else sf
    eff = c.sum(axis=0).mean() * sf  # effective library sizes
    return np.log2((c + prior) / (eff + 2 * prior) * 1e6)


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = ok.sum()
    if n == 0:
        return out
    order = np.argsort(p[ok])
    ranked = p[ok][order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.minimum(ranked, 1.0)
    out[ok] = q
    return out


def _trigamma_inverse(x: float) -> float:
    y = 0.5 + 1.0 / x
    for _ in range(50):
        tri = special.polygamma(1, y)
        d = tri * (1 - tri / x) / special.polygamma(2, y)
        y += d
        if -d / y < 1e-8:
            break
    return float(y)


def squeeze_var(s2: np.ndarray, df: float) -> tuple[np.ndarray, float, float]:
    """Smyth (2004) empirical-Bayes variance shrinkage. Returns (s2_post, df0, s0^2)."""
    s2 = np.maximum(np.asarray(s2, float), 1e-12)
    z = np.log(s2)
    e = z - special.digamma(df / 2) + np.log(df / 2)
    emean = e.mean()
    evar = np.sum((e - emean) ** 2) / max(len(e) - 1, 1) - special.polygamma(1, df / 2)
    if evar > 0:
        df0 = 2 * _trigamma_inverse(evar)
        s02 = float(np.exp(emean + special.digamma(df0 / 2) - np.log(df0 / 2)))
        post = (df0 * s02 + df * s2) / (df0 + df)
    else:
        df0, s02 = np.inf, float(np.exp(emean))
        post = np.full_like(s2, s02)
    return post, float(df0), s02


def lm_contrasts(y: np.ndarray, design: np.ndarray, contrasts: dict[str, np.ndarray],
                 moderate: bool = True) -> dict[str, pd.DataFrame]:
    """Gene-wise OLS for y (genes x samples) with named linear contrasts."""
    y = np.asarray(y, float)
    X = np.asarray(design, float)
    n, p = X.shape
    if np.linalg.matrix_rank(X) < p:
        raise ValueError("design matrix is rank deficient")
    df = n - p
    if df < 1:
        raise ValueError("no residual degrees of freedom: need more biological replicates")
    xtx_inv = np.linalg.inv(X.T @ X)
    beta = y @ X @ xtx_inv  # genes x p
    resid = y - beta @ X.T
    s2 = (resid ** 2).sum(axis=1) / df
    if moderate and len(s2) >= 50:
        s2p, df0, _ = squeeze_var(s2, df)
        df_tot = df + df0 if np.isfinite(df0) else 1e6
    else:
        s2p, df_tot = s2, df
    out = {}
    for name, c in contrasts.items():
        c = np.asarray(c, float)
        est = beta @ c
        se = np.sqrt(s2p * float(c @ xtx_inv @ c))
        t = est / se
        pv = 2 * sps.t.sf(np.abs(t), df_tot)
        tcrit = sps.t.ppf(0.975, df_tot)
        out[name] = pd.DataFrame(dict(log2fc=est, se=se, ci_lo=est - tcrit * se,
                                      ci_hi=est + tcrit * se, t=t, pvalue=pv, fdr=bh_fdr(pv),
                                      df=df_tot))
    return out


def random_effects_meta(est: np.ndarray, se: np.ndarray) -> pd.DataFrame:
    """DerSimonian-Laird random-effects meta-analysis, row-wise.

    est, se: units x studies (NaN allowed).  Returns pooled estimate, its SE, tau
    (between-study SD), a *prediction* SD for a new study/context
    (sqrt(se_pooled^2 + tau^2)), Cochran Q p-value and I^2.
    With k < 3 studies tau is very poorly estimated; ``k`` is returned so callers can flag it.
    """
    est = np.asarray(est, float)
    se = np.asarray(se, float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    w = np.where(ok, 1.0 / np.where(ok, se, 1.0) ** 2, 0.0)
    e0 = np.where(ok, est, 0.0)
    k = ok.sum(axis=1)
    sw = w.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        fixed = (w * e0).sum(axis=1) / sw
        q = (w * (e0 - fixed[:, None]) ** 2).sum(axis=1)
        c = sw - (w ** 2).sum(axis=1) / sw
        tau2 = np.where(k > 1, np.maximum((q - (k - 1)) / c, 0.0), 0.0)
        wr = np.where(ok, 1.0 / (np.where(ok, se, 1.0) ** 2 + tau2[:, None]), 0.0)
        pooled = (wr * e0).sum(axis=1) / wr.sum(axis=1)
        se_p = np.sqrt(1.0 / wr.sum(axis=1))
        i2 = np.where(q > 0, np.maximum((q - (k - 1)) / q, 0.0), 0.0)
        qp = np.where(k > 1, sps.chi2.sf(q, np.maximum(k - 1, 1)), np.nan)
    z = pooled / se_p
    pv = 2 * sps.norm.sf(np.abs(z))
    return pd.DataFrame(dict(estimate=pooled, se=se_p, tau=np.sqrt(tau2),
                             pred_sd=np.sqrt(se_p ** 2 + tau2), k=k, Q_pvalue=qp, I2=i2,
                             pvalue=pv, fdr=bh_fdr(pv)))
