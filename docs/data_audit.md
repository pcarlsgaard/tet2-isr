# Data audit and risk memo

Audit date 2026-11-12. For each series the GEO series record and every sample record were read and saved
(`results/audit/`, checksums in `results/audit/geo_reachability.tsv`); processed files were downloaded with SHA-256
(`data/raw/download_log.json`). Full fields: `data/dataset_manifest.tsv`.

## Does a direct WT/TET2 × ISR translation dataset exist?
No. A GEO DataSets query for TET2 with ribosome profiling, Ribo-seq, polysome or translation efficiency returned
0 series; TET2 with ISR/ATF4/eIF2/HRI/ISRIB returned only an unrelated valine-restriction study. Consequence:
`I_TE` cannot be estimated from public data.

## Corrections to the starting list

| Series | Brief said | GEO record says |
|---|---|---|
| GSE268127 | possible ISR translation data | RNA-seq only (ATF4 uORF1-mutant mESC, 3 vs 3) |
| GSE227766 | cycling/quiescent/senescent × thapsigargin | thapsigargin (25 nM, 3 h) only in cycling cells; Tg libraries deposited 18 months after controls, so batch is confounded with treatment. State dependence of the ISR is not estimable from it |
| GSE132090 | check replication | no biological replication: one pool of two mice per genotype, two technical libraries. Excluded |
| GSE285379 | primary human scRNA set | one library per genotype × treatment. Excluded (cells are not replicates) |
| GSE287524 | secondary | usable: CD34+ HSPC, CRISPR TET2-mutant vs WT, 3 vs 3 untreated (and 3 vs 3 LPS) |
| GSE183675 | RNA + epigenome | RNA deposited only as normalised values; LSK, 2 vs 2 per age |
| GSE65778 | (added) | confirmed: HEK293T, tunicamycin 1 h ± ISRIB, paired Ribo/RNA, 2 per arm |

## Sources used

TET2 basal contrasts (one "target" population each; they are not the same population):

| Series | Species | Target population | n WT / TET2 | Main caveat |
|---|---|---|---|---|
| GSE232384 | mouse | sorted LT-HSC (200 cells) | 3 / 4 animals | low-input noise |
| GSE212534 | mouse | sorted donor HSC, 13 months post-transplant | 2 / 2 libraries | transplant, clonal selection |
| GSE183675 | mouse | LSK, young | 2 / 2 | HSPC not HSC; normalised input |
| GSE287524 | human | CD34+ from long-term culture | 3 / 3 | cultured; edited pool |

Tet2/TET2 mRNA is not lower in the mutant samples of any study, so genotype cannot be confirmed from expression.

ISR translation (paired RNA + Ribo, ≥2 replicates per arm, found by GEO search plus the brief):

| Series | Cells | Stressor | Role in this run |
|---|---|---|---|
| GSE65778 | HEK293T | tunicamycin 1 h; + ISRIB | engaged; ISRIB arm = specificity control |
| GSE254631 / GSE254349 | U-2OS (G3BP1 rescue) | arsenite 200 µM 2 h | engaged |
| GSE282509 | HeLa | tunicamycin 16 h | engaged (late) |
| GSE227766 | IMR90 | thapsigargin 3 h | gate not passed (mean 0.57, p = 0.20) |
| GSE291029 | RPE-1 TP53-KO | WEE1 inhibitor 24 h | gate not passed |
| GSE193945 | mouse SCC | bortezomib 4.5 h; eIF2α-S51A | gate not passed; S51A = control |
| GSE129757 | LN308 | tunicamycin 2 h, 6 h | gate not passed (no reference response) |
| GSE119365 | mouse erythroblasts | iron deficiency (HRI); Hri-KO | excluded by depth rule (one library 0.14 M reads) |
| GSE256237 | HCT116 | arsenite 1 h | excluded by depth rule |

GSE119365 is the only hematopoietic and only HRI-specific translatome. It is too shallow for the rule set here;
it is the first candidate for re-quantification from raw reads.

Not used: GSE273599/600/601 (RNA-only Dmr-PERK time and dose series; needed for the late transcriptional phase,
which is not modelled), GSE307624 (ATAC tracks), GSE174764 (unreplicated doses), GSE248724/5 (basal eIF2α-S51A).

## Risks
| Risk | Consequence | Handling |
|---|---|---|
| No TET2 × ISR translatome | `I_TE` unidentified | not estimated; M1 states it as an assumption |
| Only 3 ISR-engaged contexts, all non-hematopoietic cell lines | genome-wide transfer unsupported; M2 unfit | engagement gate; M2 gated off; predictions limited to a core set |
| Reference set both gates contexts and is predicted | partly circular for those genes | flagged per row; ISRIB reversal is independent evidence |
| TET2 studies disagree | no consensus basal state | leave-heaviest-study-out replication rule |
| Small n (2–4), mixed platforms | unstable variance | moderated statistics; one column = one biological replicate |
| Legacy annotations, cross-species mapping | gene loss | one-to-one orthologs only (16,307 kept; 5,745 ambiguous pairs dropped) |
| Licences | reuse terms not individually checked | nothing redistributed; `make public` re-downloads |
| Processed counts only | no footprint-length, periodicity or rRNA QC | stated in QC table |
