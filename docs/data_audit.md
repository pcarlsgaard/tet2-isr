# Data audit and risk memo (Milestone 1 — incomplete)

## Blocker
The build container's network policy refused every primary data host: `www.ncbi.nlm.nih.gov`,
`ftp.ncbi.nlm.nih.gov`, `eutils.ncbi.nlm.nih.gov`, `zenodo.org`, `ftp.ensembl.org`, `www.ebi.ac.uk`,
`europepmc.org`, `www.nature.com`, `data.broadinstitute.org`. Reachable: PyPI, conda-forge,
`raw.githubusercontent.com`, `github.com` pages, and a web-search tool returning third-party snippets.

**Impact on inference:** no GEO sample sheet, supplementary file, licence statement or count matrix was inspected.
Fields in `data/dataset_manifest.tsv` marked `unverified` or `claimed_unverified` repeat the project brief and
have not been independently confirmed. No real-data analysis was run. Nothing in this repository is a finding
about TET2.

**To clear it:** allow the hosts above, then `make audit-data` (saves each GEO record to `results/audit/` and
counts samples) and `make download-processed`. Fill the unverified manifest fields by hand from the saved records.

## What was confirmed, and how

| Source | Confirmed | Route |
|---|---|---|
| RiboNN (`Sanofi-Public/RiboNN`) | Repo exists; human + mouse multi-task models; weights at `zenodo.org/records/17258709`; training TE comes from the paper's supplementary tables; the repo's own training file contains **fake** TEs; separate model-weights licence file exists (terms not read); predicts basal TE from 5′UTR/CDS/3′UTR sequence | README and repo page read directly |
| `CenikLab/TE_model` (Liu et al.) | Repo exists, GPL-3.0; ribo files are *not* in the repo (email request or regenerate); HeLa example on Zenodo 10594392; location of processed TE matrices not found | repo page read directly |
| GSE227766 | Title "Senescence-specific translation dysregulation desensitizes cells to stress by inhibiting activation of the integrated stress response"; human; abstract reports high eIF2α-P with suppressed ATF4 translation in senescent cells; linked article PMC11231295 | search snippet (OmicsDI). Thapsigargin/quiescent arms, RNA+Ribo pairing and replicate counts **not** confirmed |
| GSE268127 | Not confirmed. Sibling records GSE268130/GSE268132 belong to "Plasticity of the Mammalian Integrated Stress Response" (split ISR, PCK2) | search snippet |
| GSE65778 (added lead) | Title "The small molecule ISRIB reverses the effects of eIF2 phosphorylation on translation and stress granule assembly" | search snippet; design recalled from the publication, unverified |
| All TET2 series, GSE273599–601 | Nothing beyond the brief | — |

## Does a direct WT/TET2 × ISR translation-efficiency dataset exist?
None is known among the listed sources, and none was found by search. This is a negative from an incomplete
search (GEO itself was not queryable), not a proof of absence. Consequence: **`I_TE` cannot be estimated from
public data.** It can only be (a) assumed zero (M1), (b) predicted through learned state-dependence with transfer
assumptions (M2), or (c) measured.

## Provisional priority (to be revised after verification)

ISR translation (need ≥1 adequate paired RNA/Ribo contrast):
1. **GSE65778** — ISRIB arm separates eIF2α-P-dependent translation from stressor pleiotropy. Cell line.
2. **GSE227766** — paired design per brief; gives state dependence (cycling/quiescent/senescent) within one cell
   type, directly relevant to M2. Thapsigargin is pleiotropic; fibroblasts are far from HSPCs.
3. **GSE268127** — only if Ribo-seq is actually present.
Not usable for TE: GSE273599/600/601 (RNA only; use for late ATF4-dependent transcription).

TET2 basal contrasts (need ≥2 independent studies):
1. **GSE285379** (human scRNA; untreated arms only for the basal contrast; donor pseudobulk).
2. **GSE212534** (mouse HSC/MPP3 sc; animal pseudobulk; transplant caveat).
3. **GSE183675**, **GSE132090** (replication structure must be checked — pooled animals may leave no residual df,
   in which case `genotype_de` refuses to run), **GSE232384**.
GSE287524 and GSE307624 are the same study as GSE285379 and do not count as independent replication.

## Risks
| Risk | Consequence | Handling in code |
|---|---|---|
| No direct TET2 × ISR translatome | `I_TE` unidentified | M1 labels it an assumption; M2 output carries "not validated" + extrapolation flag |
| ISR translation data are from cell lines/fibroblasts | transfer to HSC is an extrapolation | held-out-cell-type split; `extrapolation_flag` when HSC state lies outside training range |
| Stressors are pleiotropic (Tg, arsenite) | signature mixes ISR and non-ISR effects | stressor class tracked; canonical vs pleiotropic reported separately; ISRIB-reversible component preferred |
| Dose not commensurable across studies | dose ↔ responsiveness confound (see identifiability memo) | both `estimated` and `metadata` dose runs; sensitivity reported |
| Pooled/technical replicates in TET2 sets | false precision | one column = one independent replicate; no df → error |
| Cell-composition shifts in TET2 samples | per-cell vs composition effects conflated | compartment-wise pseudobulk; composition exported separately |
| Licences unread | reuse terms unknown | manifest `license` column says so; no data redistributed in this repo |
