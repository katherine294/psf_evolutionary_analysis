#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --ntasks=10

module purge
module load bear-apps/2021b
module load BLAST+/2.12.0-gompi-2021b
module load BEDTools/2.30.0-GCC-11.2.0

set -euo pipefail

# FIRST -- Download ISfinder database: https://github.com/thanhleviet/ISfinder-sequences

# Define paths to directories
GENOME="/rds/homes/k/kgh742/psf_wgs_project/03.ReferenceGenomes/Psv_NCPPB3335/PsvNCPPB3335.fna"
BLAST_DB="/rds/homes/k/kgh742/psf_wgs_project/03.ReferenceGenomes/ISfinder_blast/ISfinder_blastdb"
OUTDIR="/rds/homes/k/kgh742/psf_wgs_project/03.ReferenceGenomes/Psv_NCPPB3335/ISblast_output_test"

mkdir -p "$OUTDIR"

# Define BLAST output files
BLASTOUT_FILE="$OUTDIR/allcontigs_blastout.tsv"
FILTERED_BLASTOUT="$OUTDIR/allcontigs_blastout.filtered.tsv"
FILTERED_BEDFILE="$OUTDIR/IS_hits_unmerged.bed"
FILTERED_MERGED_BEDFILE="$OUTDIR/IS_hits_merged.bed"
FINAL_BEDFILE="$OUTDIR/IS_hits_final.bed"

STRAIN="PsvNCPPB3335"

# Run BLAST
blastn \
  -query "$GENOME" \
  -db "$BLAST_DB" \
  -out "$BLASTOUT_FILE" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen" \
  -evalue 1e-10 \
  -num_threads 8

# Filter BLAST output by query coverage and percent identity
awk 'BEGIN {OFS="\t"}
{
    alignment_length = $4
    subject_coverage = (alignment_length / $13) * 100

    if (subject_coverage >= 90 && $3 >= 90) {
        print $0, subject_coverage
    }
}' "$BLASTOUT_FILE" > "$FILTERED_BLASTOUT"


# Convert to bed format
awk 'BEGIN{OFS="\t"}
{
  chrom=$1
  start=$7-1
  end=$8
  name=$2
  score=$12
  strand=($9<=$10?"+":"-")
  print chrom, start, end, name, score, strand
}' "$FILTERED_BLASTOUT" > "$FILTERED_BEDFILE"

# Merge bed file using bedtools
bedtools sort -i "$FILTERED_BEDFILE" |
bedtools merge \
  -d 20 \
  -c 4,5,6 \
  -o distinct,max,distinct \
> "${FILTERED_MERGED_BEDFILE}"

# Keep one IS family name per hit
awk 'BEGIN{OFS="\t"}{split($4,a,","); print $1,$2,$3,a[1],$5,$6}' "${FILTERED_MERGED_BEDFILE}" > "${FINAL_BEDFILE}"

# Total IS count
total_count=$(wc -l < "${FINAL_BEDFILE}")
# IS per contig
contig_count=$(cut -f1 "${FINAL_BEDFILE}"| sort | uniq -c)
# IS family counts
family_count=$(cut -f4 "${FINAL_BEDFILE}" | sort | uniq -c | sort -nr)

# Write to summary file
{
  echo "Summary of strain: ${STRAIN}"
  echo ""
  echo "Total IS count: $total_count"
  echo ""
  echo "IS count per contig:"
  echo "$contig_count"
  echo ""
  echo "IS count per family:"
  echo "$family_count"
} > "$OUTDIR/summary.txt"

