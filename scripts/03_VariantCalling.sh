#!/bin/bash
#SBATCH --job-name=snippy_align
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --output=slurm_logs/03_snippy_align_%j.out

set -euo pipefail

module purge
module load bear-apps/2019b
module load snippy/4.6.0-foss-2019b-Perl-5.30.0
module load Python/3.7.4-GCCcore-8.3.0
module load seqtk/1.3-GCC-8.3.0

PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
TRIMMED_DATA="$PROJECT_BASE/02.TrimmedReads"
SNIPPY_OUTPUT="$PROJECT_BASE/snippy_output"
REF_DIR="$PROJECT_BASE/reference_sequences"

mkdir -p "$SNIPPY_OUTPUT" slurm_logs

#### YOU REQUIRE A COMPLETE, HIGH QUALITY REFERENCE GENOME. COPY THE ASSEMBLY FASTA INTO A REFERENCE DIRCETORY #####

# Ensure chromosomal ref exists (first sequence is chromosome)
if [ ! -f "$REF_DIR/246539E_W163A3B1_chrom.fasta" ]; then
echo "Preparing chromosomal reference from hybrid fasta"
# Use seqtk to select the chromosomal contig only and save as a new file
seqkit seq --head 1 "$REF_DIR/246539E_W163A3B1.fasta" > "$REF_DIR/246539E_W163A3B1_chrom.fasta"
fi

while read -r strain; do
  [ -z "$strain" ] && continue  # skip empty lines

  echo "Processing strain: $strain"
  OUTDIR="$SNIPPY_OUTPUT/$strain"
  R1="$TRIMMED_DATA/${strain}_1P.trim.fastq.gz"
  R2="$TRIMMED_DATA/${strain}_2P.trim.fastq.gz"
  REFERENCE="$REF_DIR/246539E_W163A3B1_chrom.fasta"

  if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
    echo "Missing trimmed reads for $strain - skipping" >&2
    continue
  fi

  snippy --outdir "$OUTDIR" \
    --ref "${REFERENCE}" \
    --R1 "$R1" --R2 "$R2" \
    --mincov 10 --minfrac 0.9 --mapqual 60 \
    --unmapped --rgid "$strain" \
    --cpus ${SLURM_CPUS_PER_TASK:-4}

done < Strains.txt
