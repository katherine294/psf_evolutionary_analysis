#!/bin/bash
#SBATCH --job-name=mosdepth_summary
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/mosdepth_summary_%j.out

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2022b
module load Miniforge3/24.1.2-0

eval "$(${EBROOTMINIFORGE3}/bin/conda shell.bash hook)"
source "${EBROOTMINIFORGE3}/etc/profile.d/mamba.sh"

CONDA_ENV_PATH="/rds/homes/k/kgh742/psf_wgs_project/conda/kgh742_mosdepth"

mamba activate "$CONDA_ENV_PATH"

cd /rds/homes/k/kgh742/psf_wgs_project/Phylogenetic_analysis/MOSDEPTH_RESULTS

python /rds/homes/k/kgh742/psf_wgs_project/summarise_mosdepth.py
