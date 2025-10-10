#!/bin/bash
#SBATCH --job-name=snippy_align
#SBATCH --time=24:00:00
#SBATCH --ntasks=12
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --output=slurm_logs/03_snippy_align_%j.out

set -euo pipefail

module purge
module load bear-apps/2022a
module load Gubbins/3.3.1-foss-2022a

PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
SNIPPY_OUTPUT="$PROJECT_BASE/snippy_output"
RECOMB_DIR="$PROJECT_BASE/phylogenetic_analysis/recombination"

mkdir -p "$RECOMB_DIR"

# copy the cleaned alignment into the recombination dir (makes gubbins outputs tidy)
ALN_CLEAN="${SNIPPY_OUTPUT}/W163a3b1_core_genome_unmasked_clean.full.aln"
if [ ! -f "$ALN_CLEAN" ]; then
  echo "ERROR: cleaned alignment not found: $ALN_CLEAN"
  exit 1
fi
cp "$ALN_CLEAN" "$RECOMB_DIR/"
cd "$RECOMB_DIR"

# run Gubbins on the copied file; use basename (file is now in this directory)
ALN_BASENAME=$(basename "$ALN_CLEAN")
run_gubbins.py --prefix W163a3b1_core_gubbins "$ALN_BASENAME"

# convert Gubbins recombination GFF to BED (0-based start)
GFF="W163a3b1_core_gubbins.recombination_predictions.gff"
if [ ! -f "$GFF" ]; then
  echo "ERROR: expected Gubbins GFF not found: $GFF"
  exit 1
fi

awk '$3 == "recombination" {print $1"\t"($4-1)"\t"$5"\t"$3}' "$GFF" > recomb_regions.bed

echo "Recombination BED written to: $RECOMB_DIR/recomb_regions.bed"

