#!/bin/bash
#SBATCH --job-name=snippycore_masked
#SBATCH --time=24:00:00 # RAxML-NG can take a few horus depending on number of genomes
#SBATCH --ntasks=14
#SBATCH --nodes=1
#SBATCH --output=slurm_logs/snippycore_masked_%j.out

set -euo pipefail

# Load dependencies
module purge
module load bear-apps/2019b
module load snippy/4.6.0-foss-2019b-Perl-5.30.0
module load RAxML-NG/1.0.1-gompi-2019b

# Define directories
PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
REF_DIR="$PROJECT_BASE/reference_sequences"
SNIPPY_OUTPUT="$PROJECT_BASE/snippy_output"
MASKED_DIR="$PROJECT_BASE/Phylogenetic_analysis/masked_regions"
FINAL_SNIPPY="$PROJECT_BASE/Phylogenetic_analysis/Core_masked_alignment"

ALN_PREFIX="W163a3b1_core_genome_masked"
STRAINS=$(cat "$PROJECT_BASE/Names.txt")

mkdir -p "$FINAL_SNIPPY"
cd "$SNIPPY_OUTPUT"

# -----------------------------------------------------------------------------------------#
#                          Run Snippy-core with mask file                                  #
# -----------------------------------------------------------------------------------------#

echo "Running snippy-core with mask file..."

snippy-core \
  --prefix "$ALN_PREFIX" \
  --ref "$REF_DIR/246539E_W163A3B1_chrom.fasta" \
  --mask "$MASKED_DIR/mge_recomb_merged.bed" \
  $STRAINS

# Move output files to the final directory
mv "$SNIPPY_OUTPUT"/${ALN_PREFIX}* "$FINAL_SNIPPY"/
echo "Masked core genome alignment saved in $FINAL_SNIPPY"

# -----------------------------------------------------------------------------------------#
#                           Generate tree using RAxML-NG                                   # 
# -----------------------------------------------------------------------------------------#

cd "$FINAL_SNIPPY"

echo "Starting RAxML-NG phylogeny..."

mpirun raxml-ng-mpi \
  --all \
  --msa "${ALN_PREFIX}.aln" \
  --model GTR+G \
  --prefix "${ALN_PREFIX}_raxml" \
  --seed 2 \
  --threads 12 \
  --bs-metric fbp,tbe

echo "RAxML-NG tree construction complete: ${ALN_PREFIX}_raxml.*"












