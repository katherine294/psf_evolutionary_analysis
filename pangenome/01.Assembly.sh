#!/bin/bash
#SBATCH --job-name=spades
#SBATCH --array=1-124
#SBATCH --ntasks=4
#SBATCH --time=24:00:00
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/spades_%A_%a.out
#SBATCH --error=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/spades_%A_%a.err

set -euo pipefail

module purge
module load bear-apps/2021b
module load SPAdes/3.15.3-GCC-11.2.0

PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
TRIMMED_READS="$PROJECT_BASE/02.TrimmedReads"
OUTPATH="$PROJECT_BASE/Assemblies"
STRAIN_LIST="$PROJECT_BASE/Names.txt"

mkdir -p "$OUTPATH"

# Get strain for this array task
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$STRAIN_LIST")

R1="$TRIMMED_READS/${NAME}_1P.trim.fastq.gz"
R2="$TRIMMED_READS/${NAME}_2P.trim.fastq.gz"

# Check if both reads exist
if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "Skipping $NAME: missing read files"
    exit 0
fi

# Create output directory per strain
STRAIN_OUT="$OUTPATH/$NAME"
mkdir -p "$STRAIN_OUT"

# Run SPAdes
echo "Running SPAdes for $NAME..."
spades.py --careful \
  -1 "$R1" \
  -2 "$R2" \
  -o "$STRAIN_OUT"

echo "SPAdes completed for $NAME"
