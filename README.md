# tet2-isr: transferable translatomic rules for genotype-dependent ISR responses in TET2-deficient HSPCs

Research software and hypothesis-generation project. It asks whether transcripts and programs respond differently
to integrated stress response (ISR) activation in TET2-deficient versus WT HSPCs, and separates two mechanisms:

1. **substrate/compositional**: TET2 loss changes which mRNAs are present; a *common* ISR translation response then
   yields different output (model M1);
2. **true genotype × ISR translation interaction** `I_TE`: the stress-induced TE change itself differs (needs M2/M3
   or a direct 2×2 translatome).

## Status

| | Status |
|---|---|
| Public-data audit | Done for 24 GEO series (sample records saved in `results/audit/`); see `docs/data_audit.md` |
| Real-data run (`make public`, about 1 min after download) | Done: 4 TET2 studies, 9 paired RNA/Ribo ISR studies ingested; report in `reports/research_report.html` |
| M0 (basal TET2 contrasts), M1 (substrate-only overlay) | Run on real data |
| M2/M3 (state-dependent / low-rank) | Implemented and tested on synthetic data; **not supported by the public data** (3 ISR-engaged contexts) and not used for genotype prediction |
| Direct WT/TET2 × ISR translatome | **Does not exist in GEO** (0 series). `I_TE` is not estimated |
| Synthetic ground-truth demo (`make smoke-test`) | Kept as a method check: `reports/research_report_synthetic_demo.html`, gene IDs `SYN…` |

## What the public data show

- **No reproducible basal TET2 signature across studies.** Gene-level log2(TET2/WT) correlations between the four
  studies are −0.06 to 0.09. 902 genes reach meta-analysis FDR < 0.05, but 644 of those are carried >80% by one
  n = 2 vs 2 study; 78 survive a replication rule. No program differs consistently.
- **The canonical ISR translation response is real and specific where the ISR is engaged.** Reference uORF targets
  (ATF4, ATF5, DDIT3, IFRD1, SLC35A4 …) gain 1–3.5 log2 in TE under tunicamycin 1 h, arsenite 2 h and tunicamycin
  16 h; ISRIB removes it (set mean 2.0 → 0.4).
- **Five of eight replicated stress contexts fail that positive control** (thapsigargin 3 h borderline; bortezomib,
  WEE1 inhibitor and LN308 tunicamycin flat).
- **Genome-wide dTE does not transfer between studies.** Held-out-study Pearson is 0.01–0.03 for a pooled gene
  signature, −0.08 for sequence features, 0.18 for the nearest single context. What transfers is a small core and
  program-level structure (program recovery r ≈ 0.7–0.9).
- **Substrate-only genotype differences are negligible.** Program-level M1 interactions are all |log2| < 0.1 and
  change sign between TET2 studies; 1 of 72 intervals excludes zero.
- **Nothing in public data speaks to a genotype-specific translation effect.** That needs the experiment in
  `docs/experimental_validation_plan.md`.

## Quick start

```bash
make setup            # venv + locked dependencies
make smoke-test       # tests + synthetic end-to-end demo (no downloads)
make audit-data       # validate manifest, save GEO sample records
make public           # download processed files (1.4 GB incl. reference), build, analyse, report
```
Stages can be run singly with `CONFIG=configs/public.yml` (`analyze-isr`, `analyze-tet2`, `fit-baselines`,
`evaluate`, `report`).

## Layout

```
configs/      smoke.yml, public.yml, public_sources.tsv (files to download), prior_hypotheses.tsv
data/         dataset_manifest.tsv; raw/ and processed/ are gitignored (re-created by make public)
docs/         data_audit, assumptions_and_limits (identifiability memo), methods, model_card, decisions,
              experimental_validation_plan
src/tet2isr/  contracts, effects, stats, data_ingest/ (synthetic, geo, public adapters), tet2_transcriptome/,
              isr_translation/, modeling/, validation/, reporting/, pipeline, cli
tests/        22 tests
results/      public/ (real run), smoke/ (synthetic), audit/ (GEO records), prediction_priority_table.csv
reports/      research_report.html (real data), research_report_synthetic_demo.html, executive_summary.md
```
