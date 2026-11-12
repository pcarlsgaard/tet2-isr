# Identifiability memo, assumptions and limits

## 1. What each data type can estimate

| Quantity | Needed data | Available publicly? | Status |
|---|---|---|---|
| Basal genotype contrast in RNA, `log2 M[TET2]/M[WT]` | WT vs TET2 RNA, matched compartment, biological replicates | Yes: 4 usable studies | **Estimable** (M0), but studies disagree |
| ISR RNA response | stress vs control RNA | Yes | Estimable, per context |
| ISR `dTE` in some cell type | paired RNA + Ribo, stress vs control, replicates | Yes: 9 studies ingested, 3 contexts ISR-engaged | Estimable **in those cells**, not in HSCs |
| ISR `dTE` in WT HSCs | Ribo-seq in HSCs under ISR | Not found | Transfer assumption |
| `I_TE = dTE[TET2] − dTE[WT]` | paired RNA + Ribo in all four genotype × ISR cells | Not found | **Unidentified.** Only an assumption (M1: 0) or a model-based extrapolation (M2) |
| Linear-scale output difference | `M` per genotype + a `dTE` | via transfer | Computable under M1 assumptions; arbitrary units |
| Protein abundance, flux, survival, venetoclax sensitivity | proteomics / functional assays / a trained phenotype map | No | Out of scope; hypotheses only |

## 2. Substrate effect versus TE interaction (worked example)

`S = M × TE`. Let ISR multiply TE by the same factor `f_g` in both genotypes.

- Gene level: `log2 S_high − log2 S_basal = log2 f_g` in both genotypes → `I_TE = 0` and log-scale `I_S = 0`.
- Linear scale: `ΔS = M (f_g − 1)`, so `I_S_linear = (M_TET2 − M_WT)(f_g − 1) ≠ 0` whenever basal mRNA differs.
- **Aggregates are different.** For a pathway `P`, `log2(Σ_P M_g f_g / Σ_P M_g)` is an abundance-weighted mean, so
  it differs by genotype whenever TET2 loss changes the *relative* abundances within `P`. A non-zero pathway-level
  log interaction, a non-zero "pathway TE" interaction, or a gene-level one summed over isoforms with different
  5′UTRs can all arise with no change in any per-transcript translation rule. Likewise for shares of total
  synthesis.

`results/smoke/example_2x2_summary.csv`:

| unit | dTE WT | dTE TET2 | I_TE | ΔS WT | ΔS TET2 | I_S linear | I_S log2 |
|---|---|---|---|---|---|---|---|
| A (mRNA 10 → 20 with TET2 loss; TE ×2) | 1 | 1 | 0 | 10 | 20 | 10 | 0 |
| B (mRNA 40 = 40; TE ×0.5) | −1 | −1 | 0 | −20 | −20 | 0 | 0 |
| pathway {A,B} | −0.32 | 0 | 0.32 | −10 | 0 | 10 | 0.32 |

Rule used throughout: a genotype-specific *translation* claim requires `I_TE ≠ 0` at the transcript level.
Pathway-level or gene-level-over-isoforms interactions are reported as compositional unless that is shown.

## 3. Dose confound (found in the synthetic run)

ISR intensity is not comparable across stressors and studies, so it is estimated from the response of reference
uORF targets. If a cell state changes how strongly those same targets respond, the dose estimate absorbs the
state effect. In the simulation this removed the M2 signal entirely (predicted `I_TE` 0.00 vs true 0.47 at uORF
targets); supplying the true dose recovered it partially (0.27, interval spanning zero). Therefore:

- model-based `I_TE` needs an ISR-input measure that is independent of the translational readout (p-eIF2α, or
  ISRIB-reversible fraction);
- "enhanced response" and "same response at higher effective dose" cannot be separated from a single-dose design.

## 4. Assumptions

| # | Assumption | Where | If wrong |
|---|---|---|---|
| A1 | ISR `dTE` learned in non-hematopoietic cells transfers to HSCs | M1, M2 | all transferred predictions biased; tested only by held-out-cell-type proxy |
| A2 | Early phase: mRNA abundance unchanged by stress | M1 | late-phase outputs wrong; late transcription is not modelled |
| A3 | Basal TE identical across genotypes and set to 1 | M1 linear outputs | linear units mis-scaled per gene; log outputs unaffected |
| A4 | `dTE` scales linearly with one intensity axis | M1b, M2, M3 | thresholds/saturation (hypothesis H3) invisible |
| A5 | State dependence of `dTE` across cell types equals state dependence across genotypes | M2 `I_TE` | `I_TE` meaningless; not testable without direct data |
| A6 | Program scores from control RNA measure the relevant state (e.g. ternary-complex availability) | M2 | covariates are proxies; protein-level state unobserved |
| A7 | Ribosome occupancy change reflects initiation change | all TE | elongation stalling, uORF footprints misread as TE |
| A8 | One-to-one orthologs behave alike in mouse and human | cross-species | ambiguous pairs are dropped, not resolved |

## 5. What the public data added to the identifiability picture

- **A1 is weaker than assumed.** Between ISR-engaged contexts, genome-wide dTE agrees only weakly (Spearman
  0.14–0.25, negative between tunicamycin 16 h and arsenite). A pooled signature predicts a held-out study at
  r ≈ 0.02. Transfer is defensible for a core set (uORF targets up; ribosomal-protein mRNAs down under ER stress)
  and for program-level means, not gene by gene.
- **"Stress" is not "ISR".** Five of eight replicated stress contexts did not move the reference targets. Any
  signature built without a positive-control gate mixes ISR with unrelated effects.
- **The gate is partly circular.** Reference genes select contexts and are then predicted. Independent support:
  ISRIB abolishes their induction in the one study that tested it.
- **There is no consensus basal TET2 state** to feed a state-dependent model: state-program differences are
  within noise and change sign between studies.
- **M2 is not fit-able**: 3 contexts, 3 state covariates, and state scores that are confounded with library
  chemistry (poly-A, total RNA, 3′-tag).

## 6. Limits of the current implementation
- Differential TE uses log-CPM linear models with variance moderation, not a negative-binomial GLM. Confirm with
  DESeq2 `deltaTE` or anota2seq on real counts.
- Gene level only. Isoform/5′UTR-resolved quantification is supported by the contract (any stable ID) but no
  transcript-level quantifier is wired in.
- Read-level Ribo-seq QC (footprint length, periodicity, P-site offsets, CDS fraction, rRNA) is not implemented;
  the QC table says "not available" explicitly.
- RiboNN is not run (GPU recommended). Real sequence features are simple canonical-transcript summaries from Ensembl 113 (5′UTR length, GC, uAUG/uORF counts, TOP-like start, CDS and 3′UTR length); they predicted dTE in unseen genes at r = 0.24 within a context set and not at all across studies. Annotated TSSs are too imprecise for a real TOP call.
- Gene families are a symbol-root heuristic, not a paralogy database.
- Reference and ATF4-target gene lists were written from memory of the cited papers, then checked against the data (`reference_set_check.csv`), not copied from supplementary tables.
- M2 intervals are under-covering on held-out contexts in the synthetic run (95% nominal → ~89–90%); M3 worse
  (~81–83%). Treat intervals as optimistic.
- With k < 3 TET2 studies the between-study variance is essentially unestimated (flagged in the meta table).
- M4 (neural) not attempted: M2/M3 gave no held-out gain over M1 in the smoke test, so the brief's precondition is
  not met.
