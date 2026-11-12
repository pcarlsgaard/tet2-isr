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

| D11 | Depth rule (≥2 libraries per arm with ≥3×10⁵ assigned reads), fixed before modelling | removes HCT116 arsenite and the HRI erythroblast set; not relaxed to rescue the latter |
| D12 | ISR-engagement gate (reference-set mean dTE ≥ 0.5, p < 0.05) | 5 of 8 stress contexts show no reference response; ungated dose scaling was numerically unstable. Thapsigargin 3 h is borderline and stays out |
| D13 | M2/M3 genotype prediction disabled on public data | 3 engaged contexts; brief requires data support before advancing |
| D14 | `replicated` requires leave-heaviest-study-out significance | one n = 2 vs 2 study carried most meta-analysis hits |
| D15 | Target population chosen per study (LT-HSC, HSC, LSK, CD34+) | no common sorted population exists; they are not pooled as if identical |
| D16 | GSE132090 and GSE285379 excluded | no biological replication |
| D17 | Real report is `reports/research_report.html`; synthetic demo kept separately | avoid confusing simulated with observed |

## Unresolved scientific assumptions (not silently decided)
- Which TET2 allele classes may be pooled (KO, catalytic mutant, heterozygous CH, CRISPR). Default: never pooled.
- Which compartment stands for "HSC-enriched" across species and gating schemes.
- Whether thapsigargin/arsenite signatures should be restricted to their ISRIB-reversible component.
- Which gene sets best define state covariates; current choices are eIF2/eIF2B subunits, ribosomal proteins and
  Hallmark E2F targets.

## Open tasks
| Task | Done when |
|---|---|
| Re-quantify GSE119365 (HRI, erythroblasts) and GSE256237 from raw reads, with read-level Ribo QC | both pass the depth rule or are definitively dropped |
| Add more engaged ISR contexts (e.g. GSE174329, GSE181040, GSE102659, GSE68265, GSE69800) | ≥6 engaged contexts, enabling M2 |
| Confirm dTE with DESeq2 `deltaTE` or anota2seq | concordance table |
| Replace memory-curated gene lists with published tables | cited `programs.tsv` |
| Paired-donor model for GSE287524; LPS arm as a non-ISR stress comparator | report section |
| Late transcriptional phase from GSE273599–601 | time-resolved M1 |
| Isoform / 5′UTR-level quantification | transcript IDs through the contract |
| RiboNN basal TE as a prior | ablation row |
