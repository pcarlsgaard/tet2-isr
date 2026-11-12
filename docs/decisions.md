# Decisions and open tasks

| # | Decision | Reason |
|---|---|---|
| D1 | Build and validate on a synthetic ground-truth world first | All primary data hosts were blocked; brief stop-rule says continue smaller and record the blocker |
| D2 | Synthetic IDs are `SYN…`; no synthetic number on a real gene | Prevent simulated output being read as a finding |
| D3 | Real-gene priority table is prior-only with `NA` magnitudes | No real-data fit exists |
| D4 | Python-only stack; log-CPM linear models with EB moderation | R unavailable in the build container; DESeq2/anota2seq listed as confirmatory |
| D5 | scikit-learn removed; ridge solved in closed form | fewer dependencies; models are small |
| D6 | Dose handled two ways (estimated, metadata) and always reported side by side | the confound in assumptions §3 |
| D7 | Test-context dose treated as unknown in held-out evaluation | using test labels to set dose would leak |
| D8 | M4 not attempted | M2/M3 showed no held-out gain over M1; brief makes M4 conditional |
| D9 | Late (transcriptional) phase not modelled | needs GSE273599–601 time courses |
| D10 | GSE65778 added as a lead | ISRIB arm isolates eIF2α-P-dependent translation; unverified |

## Unresolved scientific assumptions (not silently decided)
- Which TET2 allele classes may be pooled (KO, catalytic mutant, heterozygous CH, CRISPR). Default: never pooled.
- Which compartment stands for "HSC-enriched" across species and gating schemes.
- Whether thapsigargin/arsenite signatures should be restricted to their ISRIB-reversible component.
- Which gene sets define state covariates in real data (`configs/public.yml` lists placeholders; the gene lists
  themselves have not been curated).

## Open tasks
| Task | Owner | Done when |
|---|---|---|
| Allow NCBI/Zenodo/Ensembl hosts; run `make audit-data`; fill unverified manifest fields | data scout | no `unverified` in rows used for analysis |
| Adapters GEO → contract for the top-ranked ISR and TET2 series | data scout + analysts | `contracts.load_*` passes on `data/processed/public` |
| Curate real program gene sets and ISR reference set with citations | translation analyst | `programs.tsv` for human and mouse |
| Run M0/M1 on public data; source-choice sensitivity | modeling lead | `results/public/` + report |
| Confirm dTE with DESeq2 `deltaTE`/anota2seq | translation analyst | concordance table |
| RiboNN basal TE as feature and as basal-TE prior for M1 linear outputs | modeling lead | ablation row |
| Read-level Ribo-seq QC if raw reads are pulled | translation analyst | QC columns filled |
| Transcript/5′UTR-isoform resolution | translation analyst | isoform-level IDs through the contract |
