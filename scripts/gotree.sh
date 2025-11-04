#!/bin/bash
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
