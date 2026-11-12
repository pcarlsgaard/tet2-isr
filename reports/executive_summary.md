# Executive summary

**Question.** Do transcripts and pathways respond differently to ISR activation in TET2-deficient versus WT
HSPCs, and is any difference due to a changed mRNA pool (substrate) or a changed translational response (`I_TE`)?

**Answer from public data.** They cannot answer it, and they give no positive indication of a difference.
There is no TET2 × ISR translatome in GEO. The two halves that do exist (basal TET2 transcriptomes; ISR
translatomes in other cells) were analysed and combined. The combination predicts essentially no
genotype difference from substrate effects, and it cannot address a true translational interaction.

## Findings

1. **Basal TET2 transcriptomes do not agree with each other.** Four usable studies (mouse LT-HSC, mouse
   post-transplant HSC, mouse LSK, human CD34+). Gene-level log2(TET2/WT) correlates −0.06 to 0.09 between
   studies. Of 902 meta-analysis hits, 644 are carried >80% by one 2-vs-2 study; 78 genes pass a replication
   rule. No program (proliferation, ribosomal proteins, eIF2 complex, ATF4 targets, BCL2 family, heme,
   interferon) differs consistently. Two listed datasets have no biological replication and were excluded.
2. **The canonical ISR translation response is strong where the ISR is engaged.** In tunicamycin 1 h (HEK293T),
   arsenite 2 h (U-2OS) and tunicamycin 16 h (HeLa), ATF4, ATF5, DDIT3, IFRD1 and SLC35A4 gain 1–3.5 log2 in TE.
   ISRIB reduces the set mean from 2.0 to 0.4.
3. **Most published stress conditions did not engage it.** Five of eight replicated contexts failed the
   positive control (thapsigargin 3 h borderline; bortezomib, WEE1 inhibitor, LN308 tunicamycin flat).
4. **A genome-wide "ISR translation signature" does not transfer between studies.** Held-out-study correlation is
   0.01–0.03 for a pooled signature and −0.08 for sequence features; no model beats predicting zero. A small core
   and program-level means do transfer (uORF targets +1.9 log2 synthesis potential; ribosomal proteins −0.6).
5. **Substrate-only genotype effects are negligible.** Applying the shared response to each genotype's mRNA
   pool gives program-level interactions below 0.1 log2 that change sign between TET2 studies (1 of 72 intervals
   excludes zero). Predicted direction is therefore the same in WT and TET2-deficient cells for every listed
   gene and program.
6. **State-dependent models are not supported.** Three engaged contexts cannot fit them; in simulation the same
   models also missed a planted genotype interaction unless ISR dose was known independently.

## What this means for HRI and venetoclax
Hypotheses only. Public data contain no evidence that TET2 loss changes ISR tone, HRI (EIF2AK1) expression
(log2FC −0.12, not significant), ATF4 targets or BCL2-family mRNA in a replicated way. PMAIP1 (NOXA) is higher in
TET2-deficient cells in 4 of 4 studies (meta log2FC 0.58) but fails the leave-one-study-out check. The only
HRI-specific hematopoietic translatome (GSE119365, erythroblasts) is too shallow as deposited. Nothing here
supports or refutes a link to BCL2 dependence; that needs BH3 profiling.

## Hypotheses carried forward (none favoured)
H0 no genotype-specific TE difference (default) · H1 enhanced · H2 blunted · H3 shifted dose threshold.

## Recommended next steps
1. **Minimal experiment:** WT vs Tet2-KO HSC-enriched cells, HRI-directed and non-HRI ISR inducers at three doses
   ± ISRIB; read p-eIF2α, ATF4 protein, OPP and *Atf4* mRNA; plot outputs against p-eIF2α.
2. **Decisive experiment, if step 1 shows a difference:** 2 × 2 paired RNA/Ribo, ≥4 biological replicates per
   arm, endpoint = mean `I_TE` of the core uORF set. Its data drop into the pipeline unchanged.
3. **Computational:** re-quantify GSE119365 from raw reads; add further engaged ISR contexts to reach the six
   needed for the state-dependent model.

Details: `reports/research_report.html`, `results/prediction_priority_table.csv`, `docs/`.
