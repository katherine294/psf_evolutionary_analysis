#!/bin/bash
#SBATCH --job-name=snippycore_masked
#SBATCH --time=24:00:00 # RAxML-NG can take a few horus depending on number of genomes
#SBATCH --ntasks=14
#SBATCH --nodes=1
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/snippycore_masked_%j.out

set -euo pipefail

# Load dependencies
module purge
module load bear-apps/2019b
module load snippy/4.6.0-foss-2019b-Perl-5.30.0
module load RAxML-NG/1.0.1-gompi-2019b

# Define directories and 
PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
REF_DIR="$PROJECT_BASE/reference_sequences"
SNIPPY_OUTPUT="$PROJECT_BASE/snippy_output"
MASKED_DIR="$PROJECT_BASE/Phylogenetic_analysis/masked_regions"
FINAL_SNIPPY="$PROJECT_BASE/Phylogenetic_analysis/Core_masked_alignment"

ALN_PREFIX="W163a3b1_core_genome_masked"
STRAINS=$(cat "$PROJECT_BASE/Names.txt")

if [ ! -s "$PROJECT_BASE/Names.txt" ]; then
    echo "ERROR: Names.txt not found or empty at $PROJECT_BASE"
    exit 1
fi

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

# -----------------------------------------------------------------------------------------#
#                           Collapse clades < 70 support                                   # 
# -----------------------------------------------------------------------------------------#

# Reload dependancies as they use a different bear version
module purge
module load bear-apps/2022b
module load Miniforge3/24.1.2-0

# Initialise conda/mamba
eval "$(${EBROOTMINIFORGE3}/bin/conda shell.bash hook)"
source "${EBROOTMINIFORGE3}/etc/profile.d/mamba.sh"

# Define paths
CONDA_ENV_PATH="/rds/homes/k/kgh742/psf_wgs_project/conda/${USER}_gotree"
export CONDA_PKGS_DIRS="/scratch/${USER}/conda_pkgs"

# Create environment (only needs to be done once)
if [ ! -d "${CONDA_ENV_PATH}" ]; then
    echo "Creating gotree conda environment..."
    mamba create --yes --prefix "${CONDA_ENV_PATH}" python=3.10 gotree
fi

# Activate environment
mamba activate "${CONDA_ENV_PATH}"

cd "$FINAL_SNIPPY"

# Input and output tree files from RAxML-NG
INPUT_TREE="${ALN_PREFIX}_raxml.raxml.bestTree"
OUTPUT_TREE="${ALN_PREFIX}_raxml_70BS_collapsed.tree"

# Collapse nodes below 70 bootstrap support
echo "Starting gotree clade collapse..."

gotree collapse support \
  --support 70 \
  --input "$INPUT_TREE" \
  --output "$OUTPUT_TREE"

echo "Collapsed tree saved as: $OUTPUT_TREE"







