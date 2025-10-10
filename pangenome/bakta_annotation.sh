#!/bin/bash
#SBATCH --job-name=bakta
#SBATCH --array=1-124
#SBATCH --ntasks=12
#SBATCH --time=24:00:00
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/bakta_%A_%a.out

set -euo pipefail

module purge
module load bear-apps/2022b
module load bakta/1.9.4-foss-2022b
module load PILER-CR/1.06-GCC-12.2.0

PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
ASSEMBLIES_DIR="$PROJECT_BASE/ "
OUTPATH="$PROJECT_BASE/Annotated_assemblies"
BAKTA_DB="/rds/homes/k/kgh742/psf_wgs_project/bakta_database_full/db" # Download the bakta database prior to running this script https://github.com/oschwengers/bakta
STRAIN_LIST="$PROJECT_BASE/Names.txt"

mkdir -p "$OUTPATH"

# Get strain for this array task
NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$STRAIN_LIST")
FASTA="${ASSEMBLIES_DIR}/${NAME}/contigs.fasta"      
OUTDIR="${OUTPATH}/${NAME}"

echo "[$(date)] Starting annotation for ${NAME}"

# Check input file
if [[ ! -s "$FASTA" ]]; then
  echo "ERROR: Missing or empty contigs file for ${NAME}: $FASTA"
  exit 1
fi

# Run Bakta
bakta --db "$BAKTA_DB" \
      --verbose \
      --compliant \
      --skip-crispr \
      --min-contig-length 500 \
      --genus Pseudomonas \
      --species savastanoi \
      --threads $((SLURM_NTASKS - 2)) \
      --prefix "$NAME" \
      --output "$OUTDIR" \
      "$FASTA"

echo "[$(date)] Annotation completed successfully for ${NAME}"


