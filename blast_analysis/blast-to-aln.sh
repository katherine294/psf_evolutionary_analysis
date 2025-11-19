#!/bin/bash
#SBATCH --time 6:00:00
#SBATCH --ntasks 16

module purge; module load bear-apps/2021b
module load BLAST+/2.12.0-gompi-2021b
module load BEDTools/2.30.0-GCC-11.2.0
module load Biopython/1.79-foss-2021b

set -e

#### MAKE SURE YOUR GENOMES ARE IN THE FASTA INTPUT DIRECTORY ####
# Assign names for the directories
FASTA_INPUT="full/path/Psf_genomes"
SCRIPTS="scripts"
TRIMMED_FASTA="Psf_trimmed_fasta"
BLAST_DB="blast_database"
BLAST_OUTPUT="blast_output"
FILTERED_BLAST="filtered_blast"
EXTRACTED_SEQ="extracted_sequences"
BED_FILE="bed_file"
EXTRACTED_FASTA_DIR="extracted_sequences"

QUERY="query.fasta" ##### file with the sequence of interest 

# Create directories if they don't exist
mkdir -p $TRIMMED_FASTA $SCRIPTS $BLAST_DB $BLAST_OUTPUT $FILTERED_BLAST $TOPHITS_OUTPUT $EXTRACTED_SEQ $BED_FILE $EXTRACTED_FASTA_DIR

# move scripts to new dir and make excecutable
mv trimcontig.py blast2bed.sh scripts
chmod +x scripts/*

### Make text file with basenames of input files ###
cd $FASTA_INPUT
for file in *.fasta; do base=$(basename $file .fasta); echo "${base}" >> Strains.txt; done
mv Psf_genomes.txt ../
cd ../

# Loop through strains listed in strains.txt
while read -r strain; do
  echo "Processing strain: $strain"

# Check if the input FASTA file exists
  if [ ! -f "$FASTA_INPUT/${strain}.fasta" ]; then ### check the ending for your strains is fasta not fna
echo "Error: FASTA file for strain $strain not found!"
continue
  fi

# Filter contigs so only >500bp remain using python script - could also use seqtk 
# Usage: trimcontig.py <fasta_file> <output_file>
scripts/trimcontig.py $FASTA_INPUT/${strain}.fasta $TRIMMED_FASTA/${strain}_trim.fasta

# Make blast database from Psf genomes
makeblastdb -in $TRIMMED_FASTA/${strain}_trim.fasta -parse_seqids -out $BLAST_DB/${strain}_nucl_blastdb -title "${strain}_blastdb" -dbtype nucl

# tblastn to identify protein seq in nucleotide genomes using protein effector query database
# Note: qlen is added to the end for calculating query coverage
tblastn -query "${QUERY}" -db $BLAST_DB/${strain}_nucl_blastdb -out $BLAST_OUTPUT/${strain}_blast.out \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
-evalue 1e-3
# E-value set based on visual inspection of unfiltered dataset
# around 11min cputime, 31244K memory used for one genome

# Filter blast output results by set number for query cover and percentage identity
awk '{
# Calculate alignment length and query coverage
alignment_length = ($8 - $7 + 1);
query_coverage = alignment_length / $13 * 100;
# Check if query coverage and percent ID meet the thresholds
    if (query_coverage >= 90 && $3 >=90) {
OFS="\t";  # Set the output field separator to tab
        print $0, query_coverage;
    }
}' $BLAST_OUTPUT/${strain}_blast.out > $FILTERED_BLAST/${strain}_filtered.out

# Extract chrom, sstart, end, strand dir and output into bed format
# Usage: ./blast2bed.sh <blastoutput.bls> <output.bed>
scripts/blast2bed.sh $FILTERED_BLAST/${strain}_filtered.out $BED_FILE/${strain}.bed

# Extract subject sequence using coordinates, specify -s so that getfasta runs strand aware
bedtools getfasta -name -s -fi $FASTA_INPUT/${strain}.fasta -bed $BED_FILE/${strain}.bed -fo $EXTRACTED_SEQ/${strain}_seq_extract.fasta

done < Strains.txt
echo "Processing complete."
