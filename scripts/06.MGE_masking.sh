#!/bin/bash
#SBATCH --job-name=mge_masking
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --output=slurm_logs/mge_masking_%j.out

set -euo pipefail

module purge
module load bear-apps/2021b
module load BLAST+/2.12.0-gompi-2021b
module load seqtk/1.3-GCC-11.2.0

# ----------------------------------------------------------------------------------------------------------------------
# Note: phage.bed and MGE_geneious.bed are expected to be created manually via PHASTEST and Geneious/Artemis annotations 
# ----------------------------------------------------------------------------------------------------------------------

### For phage sequences, upload chromosome sequence to PHASTEST (https://phastest.ca) ###
# Extract start and end positions of predicted prophage regions
# Save coordinates as a BED file, e.g.:
# phage.bed

### Manually Annotated mobile elements from GFF ###
# Load the chromosome GFF file into Geneious (or Artemis) and search for CDS annotated with keywords:
# "phage", "transposase", "tail", "head", "integrase",
# "terminase", "conjugal", "integrate", "integrative conjugative element",
# "conjugative transfer", "conjugative coupling"
# Export identified regions as BED file, e.g.:
#   MGE_geneious.bed

PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
REF_DIR="$PROJECT_BASE/reference_sequences"
ISFINDER_BLAST="$PROJECT_BASE/ISfinder_blast"
MGE_BED_DIR="$PROJECT_BASE/Phylogenetic_analysis/MGE"

mkdir -p "$ISFINDER_BLAST" "$MGE_BED_DIR" slurm_logs

# Convert chromosome to single-line if needed
seqtk seq "$REF_DIR/246539E_W163A3B1_chrom.fasta" > "$REF_DIR/246539E_W163A3B1_chrom_s.fasta"

# Make ISfinder DB - user must place IS_copy.fna in repo root or update path
if [ -f IS_copy.fna ]; then
  awk '/^>/{header=$0; count[header]++; if (count[header]>1){print header"_"count[header];} else {print header;} next} {print}' IS_copy.fna > "$ISFINDER_BLAST/IS_copy_unique.fna"
  makeblastdb -in "$ISFINDER_BLAST/IS_copy_unique.fna" -parse_seqids -dbtype nucl -out "$ISFINDER_BLAST/ISfinder_blastdb" -title "ISfinder_blastdb"
  blastn -query "$REF_DIR/246539E_W163A3B1_chrom_s.fasta" -db "$ISFINDER_BLAST/ISfinder_blastdb" -out "$MGE_BED_DIR/W163A3B1_chrom_isfinder_blast.txt" -outfmt 6
  awk '{print $1"\t"($7-1)"\t"$8"\t"$2}' "$MGE_BED_DIR/W163A3B1_chrom_isfinder_blast.txt" > "$MGE_BED_DIR/isfinder_blast.bed"
else
  echo "IS_copy.fna not found - please download IS elements to proceed" >&2
fi

# Concatenate each bed file into one file
echo "MGE masking step completed (check $MGE_BED_DIR for isfinder_blast.bed and manually created phage.bed / MGE_geneious.bed)"

cat "$MGE_BED/phage.bed" "$MGE_BED/isfinder_blast.bed" "$MGE_BED/MGE_geneious.bed" > "$MGE_BED/MGE_all.bed"
