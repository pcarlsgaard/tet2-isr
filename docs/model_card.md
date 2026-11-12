# Model card

| | M0 | M1 | M2 | M3 |
|---|---|---|---|---|
| Purpose | basal genotype RNA/program differences | substrate-only overlay of a common ISR dTE | state-conditioned dTE; model-based `I_TE` | low-rank residual programs |
| Real data: trained on | 4 TET2 studies | 3 ISR-engaged contexts (3 studies) + 4 TET2 studies | **not fitted for prediction** (too few contexts) | **not fitted** |
| `I_TE` output | none | 0 by assumption | none on real data | none |
| Compute | seconds, CPU | seconds | about 40 s (synthetic) | seconds |

## Held-out dTE prediction, real data
Held-out study, mean over the 3 engaged contexts (`results/public/heldout_summary.csv`). Each model is trained on
the other two studies.

| Model | Pearson | RMSE | program recovery r | 95% coverage |
|---|---|---|---|---|
| Permuted-label null | −0.01 | 0.94 | −0.22 | 0.81 |
| B0 zero | n/a | 0.80 | n/a | n/a |
| B2 reference gene set | 0.08 | 0.80 | 0.76 | n/a |
| B1 nearest context | 0.18 | 0.87 | 0.73 | n/a |
| M2a sequence only | −0.08 | 0.87 | −0.10 | 0.91 |
| M1 gene signature | 0.03 | 0.95 | 0.69 | 0.74 |
| M1b + dose | 0.01 | 0.95 | 0.72 | 0.77 |
| M2 + cell state | 0.18 | 0.91 | 0.91 | 0.77 |
| M3 low rank | 0.01 | 0.95 | 0.72 | 0.77 |

No model beats predicting zero on RMSE. Gene-level transfer fails; program-level means are recovered. M2's small
gain comes from two training contexts and is not evidence for state dependence. Intervals under-cover (74–77%
at nominal 95%). Negative control: the signature correlates 0.15–0.16 with thapsigargin and WEE1-inhibitor
contexts and −0.08 to 0.05 with LN308 tunicamycin and bortezomib. Within the pooled data, sequence features
predict unseen genes at r = 0.24 with or without family grouping.

## Synthetic ground-truth check (not biology)
Held-out-study Pearson: pooled signature 0.72, sequence only 0.62, nearest context 0.58, null 0.01; M2 0.72,
M3 0.66. Planted `I_TE` at uORF targets 0.47: M1 predicts 0; M2 predicts 0.00 with estimated dose and 0.27
[−0.15, 0.46] with true dose. The real data are far less transferable than this simulation assumed.

## Intended use
Ranking hypotheses and designing experiments. Not for claims about protein abundance, survival, self-renewal or
drug sensitivity.

## Known failure modes
Gene-level transfer across cell types; contexts labelled as stress that do not engage the ISR; dose and
responsiveness confounded; optimistic intervals; meta-analysis dominated by one small study; no genotype
validation of any kind.
