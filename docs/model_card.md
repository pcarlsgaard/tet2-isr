# Model card

| | M0 | M1 | M2 | M3 |
|---|---|---|---|---|
| Purpose | basal genotype RNA/program differences | substrate-only overlay of a common ISR `dTE` | state-conditioned `dTE`; model-based `I_TE` | low-rank residual programs |
| Parameters (smoke) | none | 1 pooled `dTE` + τ per gene | `u_g` per gene + 7 β + 3×7 γ | M1b + 2 loadings per gene + 2×4 map |
| Trained on (smoke) | — | 11 synthetic contexts, 6 studies | same | same |
| Trained on (real) | **nothing yet** | **nothing yet** | **nothing yet** | **nothing yet** |
| Informed by which data | TET2 RNA only | ISR RNA+Ribo (dTE), TET2 RNA (M) | + control RNA of ISR contexts (state), sequence features | same |
| `I_TE` output | none | 0 by assumption | prediction, flagged unvalidated | not produced |
| Compute | seconds, CPU | seconds | ~40 s with 40 bootstraps | seconds |

## Held-out performance — synthetic data only
Mean over held-out contexts, held-out-study split (`results/smoke/heldout_summary.csv`):

| Model | Pearson | RMSE | sign conc. | top-50 prec. | 95% coverage |
|---|---|---|---|---|---|
| Permuted-label null | 0.01 | 0.61 | 0.46 | 0.02 | 0.92 |
| B0 zero / uniform shift | n/a | 0.57 | n/a | n/a | n/a |
| B2 reference gene set | 0.54 | 0.49 | 0.59 | 0.43 | n/a |
| B1 nearest context | 0.58 | 0.52 | 0.95 | 0.51 | n/a |
| M2a sequence only | 0.62 | 0.46 | 0.92 | 0.52 | 0.94 |
| M1 gene signature | 0.72 | 0.41 | 0.99 | 0.60 | 0.93 |
| M1b + dose | 0.72 | 0.41 | 0.99 | 0.60 | 0.91 |
| M2 + cell state | 0.72 | 0.41 | 0.99 | 0.60 | 0.89 |
| M2c + stress class | 0.72 | 0.41 | 0.99 | 0.60 | 0.90 |
| M2d + epigenome | 0.72 | 0.41 | 0.99 | 0.60 | 0.90 |
| M3 low rank | 0.66 | 0.48 | 0.97 | 0.54 | 0.83 |

Held-out stressor and cell-type splits give the same ordering. Pleiotropic stressor in an unseen cell type is the
worst case (RMSE 0.52 vs 0.37–0.40 elsewhere). No added modality improved on M1. Gene-family grouping did not
change sequence-only accuracy here (0.776 both ways) because the simulated sequence model is low-dimensional; a
paralog-copy predictor reaches 0.92, which is the leakage a high-capacity model could exploit.

## Genotype generalisation — synthetic data only
Hidden direct 2×2 translatome, uORF-target set: observed `I_TE` 0.47 ± 0.05; M1 0; M2 0.00 (dose estimated),
0.27 [−0.15, 0.46] (true dose). Gene-level correlation of M2 `I_TE` with the direct measurement: 0.00 / 0.24.

## Intended use
Ranking hypotheses and designing experiments. Not for claims about protein abundance, survival, self-renewal or
drug sensitivity.

## Known failure modes
Dose–responsiveness confound; collinear state programs with few contexts; optimistic intervals; pleiotropic
stressors; transfer across cell types untested for HSCs; simulated truth is far simpler than biology, so every
number above is an upper bound on what real data will show.
