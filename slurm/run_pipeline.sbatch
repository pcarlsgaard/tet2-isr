#!/bin/bash
# Optional SLURM profile. CPU only; no GPU is needed for M0-M3.
#SBATCH --job-name=tet2isr
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=slurm/%x-%j.out
set -euo pipefail
CONFIG=${CONFIG:-configs/smoke.yml}
cd "${SLURM_SUBMIT_DIR:-.}"
[ -d .venv ] || make setup
# Data download is deliberately NOT run here: fetch on a login/transfer node if policy allows.
for stage in analyze-isr fit-baselines evaluate report; do make "$stage" CONFIG="$CONFIG"; done
