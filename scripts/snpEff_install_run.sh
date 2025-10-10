#!/bin/bash
#SBATCH --job-name=snpeff_pipeline
#SBATCH --time=12:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/snpeff_pipeline_%j.out

set -euo pipefail

# --------------------------------------------------------------------------------
# Load dependencies
# --------------------------------------------------------------------------------
module purge
module load bear-apps/2021b
module load snpEff/5.0e-GCCcore-11.2.0-Java-11
module load BCFtools/1.20-GCC-12.3.0

# --------------------------------------------------------------------------------
# Define directories and input files
# --------------------------------------------------------------------------------
PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
REF_DIR="$PROJECT_BASE/reference_sequences"
SNPEFF_DIR="$PROJECT_BASE/Phylogenetic_analysis/snpEff"
VCF_DIR="$PROJECT_BASE/Phylogenetic_analysis/Core_masked_alignment"

REF_FASTA="$REF_DIR/246539E_W163A3B1.fasta"
REF_GFF="$REF_DIR/246539E_W163A3B1.gff"

cd "$PROJECT_BASE"

mkdir -p "$SNPEFF_DIR"
cd "$SNPEFF_DIR"

# --------------------------------------------------------------------------------
# Step 1. Download and prepare snpEff - only required once
# --------------------------------------------------------------------------------
if [ ! -d "$SNPEFF_DIR/snpEff" ]; then
    echo "Downloading snpEff..."
    wget -q https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -O snpEff_latest_core.zip
    unzip -q snpEff_latest_core.zip
fi

cd "$SNPEFF_DIR/snpEff"

# Create custom database directory
mkdir -p data/Psavastanoibakta

# Copy reference GFF and FASTA into snpEff database structure
cp "$REF_GFF" "$SNPEFF_DIR/snpEff/data/Psavastanoibakta/"
cp "$REF_FASTA" "$SNPEFF_DIR/snpEff/data/Psavastanoibakta/"

cd "$SNPEFF_DIR/snpEff/data/Psavastanoibakta"

# Clean GFF file (remove FASTA section added by Prokka)
awk '/##FASTA/ {exit} {print}' "$(basename "$REF_GFF")" > genes.gff

# Rename FASTA for snpEff database
mv "$(basename "$REF_FASTA")" sequences.fa

# Remove the original unclean GFF to avoid confusion
rm -f "$(basename "$REF_GFF")"

# Add new genome entry to snpEff.config if not already present
if ! grep -q "Psavastanoibakta.genome" snpEff.config; then
    echo "Psavastanoibakta.genome : Pseudomonas savastanoi pv. fraxini bakta" >> snpEff.config
fi

# --------------------------------------------------------------------------------
# Step 2. Build snpEff database - only required once
# --------------------------------------------------------------------------------
echo "Building snpEff custom genome database..."
java -Xmx8g -jar snpEff.jar build -gff3 -v Psavastanoibakta

# --------------------------------------------------------------------------------
# Step 3. Annotate core SNP VCF files
# --------------------------------------------------------------------------------
cd "$VCF_DIR"

# Unmasked and masked VCFs
UNMASKED_VCF="chrom_core_nomask.vcf"
MASKED_VCF="chrom_core_mgerecomb.vcf"

# Output annotated VCFs
UNMASKED_ANN_VCF="coreSNP_annotated_nomask_bakta.vcf"
MASKED_ANN_VCF="coreSNP_annotated_allmask_bakta.vcf"

echo "Annotating unmasked and masked core SNP VCF files..."
java -Xmx8g -jar "$SNPEFF_DIR/snpEff/snpEff.jar" \
    -ud 0 -v -csvStats snpeff_summary_nomask.csv Psavastanoibakta "$UNMASKED_VCF" > "$UNMASKED_ANN_VCF"

java -Xmx8g -jar "$SNPEFF_DIR/snpEff/snpEff.jar" \
    -ud 0 -v -csvStats snpeff_summary_masked.csv Psavastanoibakta "$MASKED_VCF" > "$MASKED_ANN_VCF"

echo "Annotation complete."

# --------------------------------------------------------------------------------
# Step 4. Extract masked SNPs using bcftools
# --------------------------------------------------------------------------------
echo "Extracting MGE-only SNPs..."

# Compress and index VCFs
bgzip -f "$UNMASKED_ANN_VCF"
bgzip -f "$MASKED_ANN_VCF"
tabix -p vcf "${UNMASKED_ANN_VCF}.gz"
tabix -p vcf "${MASKED_ANN_VCF}.gz"

# Compare and extract unique SNPs in unmasked set but not in masked (i.e. MGE/recomb only)
bcftools isec -c all -n=1 -w1 \
    "${UNMASKED_ANN_VCF}.gz" "${MASKED_ANN_VCF}.gz" \
    -o MGE_recomb_unique_snps.vcf

echo "MGE-only SNP extraction complete."
echo "Output saved as: $VCF_DIR/MGE_recomb_unique_snps.vcf"

# --------------------------------------------------------------------------------
# Done
# --------------------------------------------------------------------------------
echo "snpEff re-annotation and filtering completed successfully."
