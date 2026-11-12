# Methods

Notation: gene `g`, genotype `k`, context `c` (study × cell type × stressor), `M` mRNA abundance (CPM after
median-of-ratios size factors), `TE` relative translation efficiency, `S = M × TE`.

## Ingestion and contract
All sources are converted to the tables in `src/tet2isr/contracts.py`. Checks: count columns equal the sample
sheet; no duplicate IDs; each ISR context has ≥2 biological replicates in both RNA and Ribo assays (an RNA-only
series fails). Mouse→human IDs use one-to-one orthologs only; ambiguous pairs are dropped and counted
(`data_ingest/orthologs.py`). Species are analysed separately and joined at the effect-size level.

## Public-data ingestion (`data_ingest/public.py`)
Gene universe: 19,354 protein-coding genes with an Ensembl-canonical transcript (release 113). Mouse genes map
through one-to-one MGI homology classes. Legacy UCSC IDs (GSE65778) are matched on version-less accession; RefSeq
rows are summed per gene. Sample-to-arm assignments are written out in the adapter and were read from the GEO
sample records. Depth rule: every arm needs ≥2 libraries with ≥3×10⁵ assigned reads, otherwise the context is
excluded from training and kept as descriptive.

## ISR engagement gate (`modeling.isr_engagement`)
A context defines the ISR signature only if the reference uORF-target set has mean dTE ≥ 0.5 log2 (relative to
the context median) with one-sided t-test p < 0.05 across reference genes. Thresholds are in `configs/public.yml`.
Contexts that fail are scored against the signature as a negative control (`unengaged_context_check.csv`). Each
stress signature is also compared with its ISR-blocked counterpart (ISRIB, eIF2α-S51A, Hri-KO;
`isr_specificity.csv`).

## Basal TET2 contrasts (`tet2_transcriptome/analysis.py`)
Single-cell data are summed to sample × compartment pseudobulks (min 10 cells). Per study and compartment:
`log2CPM ~ genotype`, OLS with empirical-Bayes variance moderation (Smyth 2004), BH FDR. One column is one
independent biological replicate; fewer than 2 per genotype is an error. Composition (compartment fractions per
sample) is exported separately. Across studies: DerSimonian–Laird random effects per gene, reporting τ, I², Q and
the number of studies agreeing in sign. Because inverse-variance weights can be dominated by one small study,
each gene also gets the heaviest study's weight share and a leave-heaviest-study-out estimate; `replicated` means
FDR < 0.05, k ≥ 3, same sign in ≥3 studies, and p < 0.05 without the heaviest study.

## Differential TE (`isr_translation/dte.py`)
Per context: `log2CPM ~ stress + assay + stress:assay`, assays normalised separately. Reported separately with SE
and 95% CI: RNA change (`stress`), ribosome-occupancy change (`stress + stress:assay`), `dTE` (`stress:assay`).
Same contrast structure as DESeq2 `deltaTE`. Direct `I_TE`, when a 2×2 translatome exists, is the three-way
`genotype:stress:assay` coefficient. `dTE` is relative (centred within context).
Output: `results/<run>/isr_translation_signatures.parquet`.

## Models (`modeling/models.py`)
- **M0** per-sample program scores (mean log2CPM of members minus sample mean), Welch test by genotype.
- **M1** pooled `dTE` per unit intensity (random effects across contexts; predictive SD includes between-context
  τ) applied to both genotypes' basal `M`. Outputs: linear `ΔS` per genotype and `I_S_linear` (arbitrary units,
  basal TE ≡ 1), program-level log interaction from abundance-weighted aggregation, share-of-total interaction.
  Monte-Carlo intervals over `dTE ~ N(mean, pred_sd)` and bootstrap of biological replicates. `I_TE ≡ 0`.
- **M2** `dTE[g,c] = I_c ( u_g + Σ_j z_cj · [1, X_g]·γ_j )`. `u_g` gene effect shrunk to a ridge sequence model
  `X_g β` (used alone for unseen genes); `z_c` standardised basal program scores from control RNA; `I_c`
  intensity. Alternating ridge fits. Genotype prediction: `I_TE[g] = I · Σ_j (z_TET2,j − z_WT,j)·[1, X_g]·γ_j`.
  Uncertainty: joint bootstrap over ISR source studies and TET2 biological replicates; sign stability; rank
  stability (mean pairwise Spearman across bootstraps); extrapolation flag when a genotype's state lies outside
  the training range.
- **M3** rank-2 SVD of the residual gene × context matrix; context scores ridge-regressed on `z_c`. Stability by
  delete-2-contexts jackknife. The `+epigenome` ablation adds `epi_*` gene features to M2.
- **Support rule.** M2/M3 genotype predictions run only with at least `min_contexts_for_state_model` engaged
  contexts (6 in the public config). The public run has 3, so they are not produced.
- **M4** not implemented (precondition unmet; `docs/decisions.md`).

Intensity: `estimated` = mean `dTE` of the reference uORF-target set minus context median, scaled to mean 1;
`metadata` = stated dose. Both are run for genotype predictions.

## Validation (`validation/evaluate.py`)
Splits (need ≥2 training studies): held-out stressor (dose variants grouped), held-out cell type, held-out study; assertions forbid shared
study / cell type / stressor between train and test. Test-context dose is treated as unknown (training mean).
Baselines: zero/uniform shift, nearest context in state space, reference gene set, permuted labels. Metrics:
Pearson, Spearman, MAE, RMSE, sign concordance among significant observed genes, top-50 precision, program-level
recovery, coverage at 50/80/95% and interval width. Held-out-stressor results are also split by canonical vs
pleiotropic and by whether the test cell type was seen. Gene-family leakage: sequence-only model with random vs
family-grouped folds, plus a paralog-copy upper bound. Genotype generalisation: predictions compared with a
reserved direct 2×2 translatome at gene and program level (synthetic only at present).

## Intersection analysis (`intersection.py`)
Pairwise partial Spearman between each TET2 contrast and each context's dTE (covariates: WT expression, CDS
length, CDS and 5′UTR GC, 3′UTR length). TET2 cores by random-effects meta-analysis with the replication rule
(independent studies) or by averaging populations within a study (shared animals; no pooling of evidence). Set
tests compare a core's mean dTE with 2,000 expression-decile-matched random sets. Overlap with the ISR translation
core by hypergeometric test. ISR-blocked and gate-failed contexts are the empirical null.

## Reproducibility
`requirements.lock`, `environment.yml`, `Dockerfile`; seeds in configs; `results/<run>/provenance.json` holds
config, software versions, git commit and SHA-256 of every input and output. Downloads log URL, time and checksum.
