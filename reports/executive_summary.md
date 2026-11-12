# Executive summary

**Question.** Do transcripts and pathways respond differently to ISR activation in TET2-deficient versus WT
HSPCs, and is any difference due to a changed mRNA pool (substrate) or to a changed translational response
(`I_TE`)?

**What exists now.** A tested, reproducible pipeline (data contract → replicate-aware TET2 contrasts →
differential TE → models M0–M3 → held-out validation → report) and a set of method findings from a synthetic
ground-truth run. **No public data have been analysed**: the build environment could not reach GEO, Zenodo,
Ensembl or journal sites, so dataset metadata are unverified and there are no real-data predictions. Every number
below is from simulation and describes the method, not TET2.

## Findings about what can and cannot be learned

1. **Public data cannot identify `I_TE`.** No direct WT/TET2 × ISR translatome was found among the listed
   sources. Separate TET2 RNA and ISR Ribo-seq datasets support only (a) basal genotype contrasts, (b) ISR
   translation responses in other cell types, and (c) a substrate-only overlay (M1).
2. **A substrate effect alone produces apparent interactions at pathway level.** With identical per-gene `dTE`,
   pathway-level log synthesis and "pathway TE" differ by genotype whenever relative abundances within the pathway
   differ (worked example: 0.32 log2 from composition alone). Gene-level-over-isoform summaries behave the same
   way. Only transcript-level `I_TE` supports a translational claim.
3. **Transfer of a pooled ISR signature works in simulation; extra modalities add nothing.** Held-out-study
   Pearson: pooled signature 0.72, sequence-only 0.62, nearest context 0.58, null 0.01. Cell-state interactions,
   stress class, an epigenome feature and a low-rank model gave no gain (0.72, 0.72, 0.72, 0.66).
4. **Good held-out `dTE` prediction says nothing about genotype interaction.** The simulation contained a true
   `I_TE` of +0.47 log2 at uORF targets. M1 predicts zero by assumption; M2 predicted 0.00 when ISR dose was
   inferred from marker genes and 0.27 (interval −0.15 to 0.46) only when the true dose was given.
5. **Dose and responsiveness are confounded.** Without an independent measure of ISR input (p-eIF2α),
   "TET2-deficient cells respond more" and "TET2-deficient cells received a larger effective dose" are
   indistinguishable.
6. **Per-gene detection is underpowered; gene-set detection is not.** Three replicates per arm found no individual
   gene at FDR 5% but detected the uORF-target set clearly.

## Hypotheses carried forward (none favoured)
H0 no genotype-specific TE difference (default) · H1 enhanced adaptive ISR · H2 blunted adaptive ISR ·
H3 shifted dose threshold. They differ in dose-response shape at matched p-eIF2α, so a dose series separates them.

Relevance to HRI and venetoclax, stated as hypotheses only: Tet2-deficient cells' reported dependence on
NCOA4-mediated ferritinophagy and mitochondrial iron makes a heme- or mitochondria-linked route to HRI
plausible, but **there is no evidence in this project that TET2 loss activates HRI**, and iron dependence may be
independent of eIF2α kinases. A link from ISR output to BCL2-family dependence would need BH3 profiling and a
venetoclax dose-response; it cannot be inferred from predicted translation.

## Recommended next steps
1. Allow the data hosts and complete the audit (`make audit-data`, `make download-processed`); write adapters for
   the top-ranked series (provisionally GSE65778 and GSE227766 for ISR translation; GSE285379 and GSE212534 for
   TET2). Then run M0/M1 on real data.
2. Run the minimal experiment: WT vs Tet2-KO HSC-enriched cells, HRI-directed and non-HRI ISR inducers at three
   doses ± ISRIB; read p-eIF2α, ATF4 protein, OPP and *Atf4* mRNA. Plot outputs against p-eIF2α.
3. If genotypes differ at matched p-eIF2α, run the 2 × 2 paired RNA/Ribo experiment with ≥4 biological replicates
   per arm and a pre-specified uORF-target-set endpoint. Its data drop into `tet2_isr_holdout/` unchanged.

Details: `docs/experimental_validation_plan.md`, `docs/assumptions_and_limits.md`, `docs/model_card.md`,
`results/prediction_priority_table.csv`.
