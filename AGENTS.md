# AGENTS.md — rules for anyone (human or agent) working in this repository

## Hard rules (from the project brief, section K)
- Human TET2 + LPS data are not TET2 + ISR translation data.
- RNA-seq-only data are never paired ribosome-profiling data (`contracts.check_paired_rna_ribo` enforces this).
- Predicted TE / ribosome occupancy is not protein abundance. ATF4 mRNA is not ATF4 translational output.
- Multiplying separate fold changes does not identify a TET2 × ISR interaction. M1 writes `I_TE = 0` with the
  note "zero by M1 assumption, not an estimate".
- Technical replicates are never biological replicates. Single cells are never replicates (pseudobulk first).
- Do not state that TET2 loss activates HRI.
- Latent factors are descriptive, not causal.
- Negative controls, held-out sets and unusable sources stay in the report.
- Blocked download or tool: continue smaller, record the blocker and its effect on inference.

## Added from the public-data run
- A stress condition is not an ISR condition until the reference targets respond (engagement gate).
- A meta-analysis hit carried by one small study is not a replicated finding.
- Do not relax the depth rule or the gate to admit a dataset you would like to use.

## Data contract
`src/tet2isr/contracts.py` is the single interface between ingestion and analysis. A new source (GEO series or
in-house WT/Tet2/HRI/VEN data) is added by writing an adapter that emits those tables; no analysis code changes.
A direct genotype × ISR translatome goes in `tet2_isr_holdout/` and is read only by the validation stage.

## Evidence labels
Every output row carries `evidence` in {observed, transferred, speculative}. Synthetic gene IDs start with `SYN`.
Never attach a synthetic number to a real gene symbol.

## Roles and ownership
| Role | Owns | Acceptance criterion |
|---|---|---|
| Evidence/data scout | `data/dataset_manifest.tsv`, `docs/data_audit.md`, `data_ingest/` | every field has a verification status; `make audit-data` passes |
| TET2 analyst | `tet2_transcriptome/` | replicate-aware contrasts per compartment; ≥2 independent studies meta-analysed with heterogeneity |
| Translation analyst | `isr_translation/`, `qc/` | RNA, Ribo and TE changes reported separately with SE; QC gaps stated |
| Modeling lead | `modeling/` | M0–M3 run on identical splits; uncertainty from study-level resampling |
| Validation reviewer | `validation/`, `tests/`, `docs/model_card.md` | leakage assertions pass; null model at chance; coverage reported |
| Integration lead | `pipeline.py`, `reporting/`, `Makefile` | `make smoke-test` regenerates every table and figure |

One owner per source file; interface changes go through `contracts.py` in a reviewed commit.

## Open tasks
See `docs/decisions.md` → "Open tasks".
