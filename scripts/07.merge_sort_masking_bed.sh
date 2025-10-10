#!/bin/bash
#SBATCH --job-name=merge_beds
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --nodes=1
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/merge_beds_%j.out

set -euo pipefail

module purge
module load bear-apps/2023a
module load BEDTools/2.31.1-GCC-12.3.0

PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
MGE_DIR="$PROJECT_BASE/Phylogenetic_analysis/MGE"
RECOMB_DIR="$PROJECT_BASE/Phylogenetic_analysis/recombination"
OUTPUT_DIR="$PROJECT_BASE/Phylogenetic_analysis/masked_regions"

mkdir -p "$OUTPUT_DIR"

# required files
MGE_BED="$MGE_DIR/MGE_all.bed"
RECOMB_BED="$RECOMB_DIR/recomb_regions.bed"
if [ ! -f "$MGE_BED" ]; then echo "ERROR: missing MGE file: $MGE_BED"; exit 1; fi
if [ ! -f "$RECOMB_BED" ]; then echo "ERROR: missing recomb file: $RECOMB_BED"; exit 1; fi

cat "$MGE_BED" "$RECOMB_BED" > "$OUTPUT_DIR/mge_recomb_combined.bed"
bedtools sort -i "$OUTPUT_DIR/mge_recomb_combined.bed" | bedtools merge -i - > "$OUTPUT_DIR/mge_recomb_merged.bed"

echo "Final masking file created: $OUTPUT_DIR/mge_recomb_merged.bed"
