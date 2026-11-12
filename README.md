# tet2-isr: transferable translatomic rules for genotype-dependent ISR responses in TET2-deficient HSPCs

Research software and hypothesis-generation project. It asks whether transcripts and programs respond differently
to integrated stress response (ISR) activation in TET2-deficient versus WT HSPCs, and separates two mechanisms:

1. **substrate/compositional** — TET2 loss changes which mRNAs are present; a *common* ISR translation response then
   yields different output (M1);
2. **true genotype × ISR translation interaction** `I_TE` — the stress-induced TE change itself differs (M2/M3, or a
   direct 2×2 translatome).

## Current status — read this first

| | Status |
|---|---|
| Pipeline (ingest contract → TET2 contrasts → differential TE → M0–M3 → held-out validation → report) | **Implemented, tested, runs end to end** on a synthetic ground-truth dataset (`make smoke-test`, ~1 min, no downloads) |
| Public-data audit | **Partial.** The build session could not reach NCBI/GEO, Zenodo, Ensembl or journal sites. Sample-level metadata are *unverified*; see `docs/data_audit.md` |
| M0/M1 on real public data | **Not run** (no data could be downloaded). `make download-processed` and `CONFIG=configs/public.yml` stages are wired but contingent on network access and per-study adapters |
| Real-gene predictions | Prior-knowledge hypothesis table only (`results/prediction_priority_table.csv`), with **no fitted magnitudes** |
| Direct WT/TET2 × ISR translatome | None known in public data. `I_TE` is unidentified from public sources |

Everything under `results/smoke/` and in `reports/research_report.html` is **simulated** (gene IDs `SYN…`). It
demonstrates the method and its limits; it is not evidence about TET2.

## Quick start

```bash
make setup            # venv + locked dependencies (requirements.lock)
make smoke-test       # tests + synthetic data -> M0..M3 -> validation -> reports/research_report.html
make audit-data       # validates the manifest; fetches GEO records if NCBI is reachable
make download-processed            # contingent: needs www.ncbi.nlm.nih.gov + ftp.ncbi.nlm.nih.gov
make analyze-isr analyze-tet2 fit-baselines evaluate report CONFIG=configs/public.yml   # contingent on data
```

Individual stages (`analyze-isr`, `fit-baselines`, `evaluate`, `report`) default to `configs/smoke.yml`.

## What the synthetic run already shows (method findings, not biology)

- **2×2 example** (`results/smoke/example_2x2_summary.csv`): shared `dTE` + different basal mRNA gives a linear-scale
  interaction of 10 units with `I_TE = 0` and gene-level log `I_S = 0`. A two-gene *pathway* nevertheless shows a
  log-scale interaction of 0.32, purely from abundance weighting.
- **Held-out prediction**: a pooled gene signature (M1) predicts `dTE` in unseen studies with mean Pearson 0.72;
  sequence-only 0.62; nearest-context 0.58; permuted-label null 0.01. Adding cell-state interactions (M2), stress
  class or an epigenome feature does not improve held-out `dTE` (0.72), and the low-rank M3 is worse (0.66).
- **Held-out `dTE` accuracy cannot detect a genotype interaction.** The simulation contains a real `I_TE` (+0.47 log2
  at uORF targets, measured in a hidden 2×2 translatome). M1 predicts 0 by assumption. M2 predicts 0.00 when ISR
  dose is estimated from marker genes and 0.27 (95% interval −0.15 to 0.46) only when the true dose is supplied.
  **ISR dose and state-dependent responsiveness are confounded unless dose is anchored independently
  (p-eIF2α).**
- **Power**: with 3 biological replicates per arm, no single gene reached FDR < 0.05 for `I_TE` in the hidden 2×2
  experiment; the uORF-target *set* was clearly detected. Gene-set-level `I_TE` is the realistic readout.

## Layout

```
configs/      smoke.yml, public.yml, prior_hypotheses.tsv
data/         dataset_manifest.tsv; raw/ and processed/ are gitignored
docs/         data_audit, identifiability memo (in assumptions_and_limits), methods, model_card, decisions,
              experimental_validation_plan
src/tet2isr/  contracts, effects, stats, data_ingest/, tet2_transcriptome/, isr_translation/, modeling/,
              validation/, reporting/, pipeline, cli
tests/        19 tests: effect arithmetic, count alignment, ortholog joins, pseudobulk, leakage, model shapes
results/      smoke/ (synthetic run + provenance.json), prediction_priority_table.csv (prior-only)
reports/      research_report.html, executive_summary.md
slurm/        optional batch profile
```

## Documents
`docs/data_audit.md` · `docs/assumptions_and_limits.md` · `docs/methods.md` · `docs/model_card.md` ·
`docs/decisions.md` · `docs/experimental_validation_plan.md` · `reports/executive_summary.md` · `AGENTS.md`
