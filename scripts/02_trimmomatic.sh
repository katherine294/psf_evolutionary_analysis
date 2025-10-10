#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --output=slurm_logs/02_trimmomatic_%j.out

set -euo pipefail

module purge
module load bear-apps/2019b
module load Java/11
module load Trimmomatic/0.39-Java-11

PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
DATADIR="$PROJECT_BASE/01.RawData"
TRIMMED_DATA="$PROJECT_BASE/02.TrimmedReads"
ADAPTOR_SEQ="$PROJECT_BASE/01.RawData/TruSeq2-PE.fa" #This file is saved in 

mkdir -p "$TRIMMED_DATA" slurm_logs

THREADS=${SLURM_CPUS_PER_TASK:-4}

cd "$DATADIR"

# If names.txt exists use it, otherwise detect pairs
if [ -f names.txt ]; then
  mapfile -t SAMPLES < names.txt
else
  SAMPLES=()
  for f in *_1.fastq.gz; do [ -e "$f" ] || continue; SAMPLES+=("$(basename "$f" _1.fastq.gz)"); done
fi

for base in "${SAMPLES[@]}"; do
  echo "Processing sample: $base"
  IN1="$DATADIR/${base}_1.fastq.gz"
  IN2="$DATADIR/${base}_2.fastq.gz"
  OUT1P="$TRIMMED_DATA/${base}_1P.trim.fastq.gz"
  OUT1U="$TRIMMED_DATA/${base}_1U.trim.fastq.gz"
  OUT2P="$TRIMMED_DATA/${base}_2P.trim.fastq.gz"
  OUT2U="$TRIMMED_DATA/${base}_2U.trim.fastq.gz"

  if [ ! -f "$IN1" ] || [ ! -f "$IN2" ]; then
    echo "Skipping $base - input fastq not found" >&2
    continue
  fi

  java -jar "$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar" PE -threads "$THREADS" \
    "$IN1" "$IN2" \
    "$OUT1P" "$OUT1U" \
    "$OUT2P" "$OUT2U" \
    ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2:True \
    LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:50

done

echo "Trimming complete. Results saved to: $TRIMMED_DATA"
