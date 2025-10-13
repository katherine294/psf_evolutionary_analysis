#!/bin/bash
#SBATCH --job-name=panaroo
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=6750M
#SBATCH --ntasks=54
#SBATCH --nodes=1
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/panaroo_%A.out
#SBATCH --error=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/panaroo_%A.err

set -euo pipefail

# Load modules
module purge
module load bear-apps/2022a
module load panaroo/1.5.0-foss-2022a

# Define paths
PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
ANNOTATION_DIR="$PROJECT_BASE/Annotated_assemblies"
OUTPATH="$PROJECT_BASE/Panaroo_output"
GFF_INPUT="$OUTPATH/Panaroo_gff_input"

mkdir -p "$OUTPATH" "$GFF_INPUT"


# Collect all Bakta .gff3 and .fna files

echo "Collecting annotated files into Panaroo input directory..."

cd "$ANNOTATION_DIR"

# Copy all relevant annotation outputs
find . -maxdepth 2 -type f -name "*.gff3" -exec cp {} "$GFF_INPUT" \;
find . -maxdepth 2 -type f -name "*.fna" -exec cp {} "$GFF_INPUT" \;

# Create a Panaroo input list file

echo "Preparing Panaroo input list..."

cd "$GFF_INPUT"

# Ensure only matching strain files are paired correctly
for gff in *.gff3; do
    base=$(basename "$gff" .gff3)
    if [[ -f "${base}.fna" ]]; then
        echo "$GFF_INPUT/${base}.gff3 $GFF_INPUT/${base}.fna"
    else
        echo "Warning: No FASTA found for ${base}" >&2
    fi
done > "$OUTPATH/Panaroo_input.txt"

# -----------------------------------------------------------------------------
# Run Panaroo
# -----------------------------------------------------------------------------
echo "Starting Panaroo analysis..."

cd "$OUTPATH"

panaroo \
    -i "$OUTPATH/Panaroo_input.txt" \
    -o "$OUTPATH/Panaroo_Psf" \
    -t "$SLURM_NTASKS" \
    --clean-mode strict \
    -c 0.98 \
    -f 0.7 \
    --len_dif_percent 0.98 \
    --aligner mafft \
    --core_threshold 0.98 \
    -a pan \
    --merge_paralogs

echo "Panaroo completed successfully at $(date)"

