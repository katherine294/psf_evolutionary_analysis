#!/bin/bash
#SBATCH --job-name=snippy_align
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/03_snippy_align_%j.out

set -euo pipefail

module purge
module load bear-apps/2019b
module load snippy/4.6.0-foss-2019b-Perl-5.30.0

PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
SNIPPY_OUTPUT="$PROJECT_BASE/snippy_output"
REF_DIR="$PROJECT_BASE/reference_sequences"

mkdir -p "$SNIPPY_OUTPUT" slurm_logs

ALN_PREFIX="W163a3b1_core_genome_unmasked"
STRAINS=$(cat Names.txt)

cd $SNIPPY_OUTPUT

snippy-core \
  --prefix "$ALN_PREFIX" \
  --ref "$REF_DIR/246539E_W163A3B1_chrom.fasta" \
  "$STRAINS"

ALN_CAND="$SNIPPY_OUTPUT/${ALN_PREFIX}.full.aln"

# Clean alignment (replace ambiguous or low-quality sites with "X")
snippy-clean_full_aln $ALN_CAND > $SNIPPY_OUTPUT/${ALN_PREFIX}_clean.full.aln
echo "Cleaned alignment written to: ${SNIPPY_OUTPUT}/${ALN_PREFIX}_clean.full.aln"

echo "Per-sample snippy complete"
