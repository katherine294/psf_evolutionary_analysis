#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --array=1-169
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/snippy_output/MOSDEPTH_RESULTS/slurmfiles/slurm-%A_%a.out

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2022b
module load Miniforge3/24.1.2-0

eval "$(${EBROOTMINIFORGE3}/bin/conda shell.bash hook)"
source "${EBROOTMINIFORGE3}/etc/profile.d/mamba.sh"

CONDA_ENV_PATH="/rds/homes/k/kgh742/psf_wgs_project/conda/kgh742_mosdepth"

if [ ! -d "$CONDA_ENV_PATH" ]; then
  echo "Creating mosdepth environment..."
  mamba create --yes --prefix "$CONDA_ENV_PATH" python=3.10 mosdepth
fi

mamba activate "$CONDA_ENV_PATH"

INPATH="/rds/homes/k/kgh742/psf_wgs_project/snippy_output"
OUTPATH="${INPATH}/MOSDEPTH_RESULTS"
STRAIN_LIST="${INPATH}/Strains.txt"

mkdir -p "$OUTPATH/slurmfiles"

NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$STRAIN_LIST")
BAM="${INPATH}/${NAME}/snps.bam"

if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found for ${NAME}"
  exit 1
fi

echo "Running mosdepth for ${NAME}"
mosdepth -t 4 -n -b 500 "$OUTPATH/${NAME}" "$BAM"

echo "Mosdepth completed for ${NAME}"
