# Experimental validation proposal (design only; nothing here has been run)

What the public-data run changes: (i) confirm ISR engagement in every sample with a positive control (ATF4 protein
or reference-set TE), because most published "stress" conditions did not engage it; (ii) pre-specify the endpoint
on the core set that transferred across studies (ATF4, ATF5, DDIT3, IFRD1, SLC35A4, PTP4A1, PPP1R15B), not on
genome-wide signatures; (iii) include an ISRIB arm, which was the cleanest specificity control available.

## Single most informative falsification experiment
**2 × 2 paired RNA-seq + Ribo-seq (or polysome-seq) in HSC-enriched cells: WT vs Tet2-KO × vehicle vs ISR
inducer, at an early time point, with ISR input matched by p-eIF2α.**
Pre-specified primary endpoint: mean `I_TE` over a fixed uORF-target set (three-way genotype × stress × assay
interaction). Outcomes: interval excludes zero → M1's "no interaction" assumption is false; interval inside a
pre-set equivalence margin (e.g. ±0.25 log2) → H0 stands and genotype differences in output are substrate-driven.
This is the only design that identifies `I_TE`; it also yields the genotype-specific `M` needed for the M1
substrate calculation, so both mechanisms are separated in one experiment.

## Minimal targeted first experiment (TET2 / HRI / ATF4 context)
Sorted LSK or HSC-enriched cells, WT vs Tet2-KO, short ex vivo culture:
- arms: vehicle; HRI-directed activator (e.g. BTdCPU) at 3 doses; a non-HRI inducer as kinase control
  (e.g. halofuginone for GCN2); each ± ISRIB;
- readouts, same cells, 1–2 h: p-eIF2α and total eIF2α (flow or capillary western), ATF4 protein, ATF4 5′UTR
  translational reporter if available, OPP incorporation, *Atf4* mRNA by qPCR;
- interpretation: plot ATF4 protein and OPP against p-eIF2α. Same curve in both genotypes → no interaction (H0).
  Vertical shift at matched p-eIF2α → altered translational responsiveness (H1/H2). Horizontal shift in
  p-eIF2α per dose → altered kinase input or feedback (H3), which is upstream of translation rules.
  *Atf4* mRNA separates substrate from TE for this one gene. ISRIB reversal confirms eIF2α-P dependence.
No sequencing needed; decides whether the 2 × 2 translatome is worth its cost and fixes its dose.

## Factorial design for the full study
| Factor | Levels |
|---|---|
| Genotype | WT, Tet2-KO (littermates; optionally Tet2+/− to model CH) |
| ISR input | vehicle, low, moderate, high (≥4 points to resolve threshold vs amplitude) |
| Time | early 1–2 h (translational), late 6–8 h (ATF4-dependent transcription) |
| Kinase axis (secondary) | HRI-directed vs GCN2- or PERK-directed inducer; Hri-KO or kinase inhibitor as specificity control |
| Venetoclax (secondary, separate arm) | ± venetoclax on top of fixed ISR dose; never as a substitute ISR inducer |

Controls: vehicle per genotype and time; ISRIB rescue; unstressed sorted-population purity and compartment
composition; cell-number-matched inputs; spike-ins if global translation changes are to be read from sequencing.
Replication: independent animals (or independent pools from distinct animals), never repeat libraries from one
sample. Randomise treatment within animal across wells; blind genotype for flow gating, western quantification
and library prep order; balance genotype across prep batches and lanes.

## Staged measurement panel
1. p-eIF2α, ATF4 protein/reporter, OPP (cheap; dose-finding; go/no-go).
2. Target protein abundance for ~10 predicted genes; RNA-seq at early and late time.
3. Targeted or genome-wide translatome (low-input Ribo-seq or polysome profiling) at the matched dose.
4. Functional assays (competitive transplantation, serial replating/self-renewal, apoptosis/BH3 profiling,
   lineage output) **only if** stages 1–3 show a genotype difference.

## Power (planning numbers, to be replaced by pilot variance)
Set-level endpoint. If the per-gene SE of `I_TE` with `n` replicates per arm is `σ₃·√(3/n)` and the set has `m`
effectively independent genes, the detectable effect at 80% power, two-sided α = 0.05, is about
`2.8 · σ₃·√(3/n) / √m`. The simulation gave `σ₃ ≈ 0.49` log2 at 3 replicates; real low-input HSC libraries are
likely noisier, and gene–gene correlation reduces `m`.

| replicates per arm | m = 10 | m = 20 | m = 50 |
|---|---|---|---|
| 3 | 0.43 | 0.31 | 0.19 |
| 4 | 0.38 | 0.27 | 0.17 |
| 6 | 0.31 | 0.22 | 0.14 |

(minimum detectable set-mean `I_TE`, log2). Per-gene discovery is not realistic at these sizes: in the
simulation, 3 replicates detected zero individual genes at FDR 5% despite a true 0.4 log2 effect in 6% of genes.
For stage 1 (flow/western), a 1.5-fold genotype difference in ATF4 at matched p-eIF2α with a CV of 25% needs
about 7 animals per genotype (two-sample t, 80% power).

## Readouts that separate substrate from TE
| Observation | Substrate effect | TE interaction |
|---|---|---|
| Basal mRNA differs; protein/RPF fold change under ISR equal | yes | no |
| mRNA equal; RPF/mRNA fold change under ISR differs at matched p-eIF2α | no | yes |
| Different 5′UTR isoform usage by genotype | yes (isoform substrate) | appears as gene-level TE change; resolve with isoform-level quantification |
| Difference disappears when plotted against p-eIF2α instead of dose | upstream signalling | no |
