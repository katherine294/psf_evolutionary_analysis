#!/bin/bash
#SBATCH --job-name=download_fastqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/01_download_fastqc_%j.out

set -euo pipefail

module purge
module load bear-apps/2019b
module load FastQC/0.11.9-Java-11
module load MultiQC/1.9-foss-2019b-Python-3.7.4

# --- Edit this base path if different ---
PROJECT_BASE="/rds/homes/k/kgh742/psf_wgs_project"
DATADIR="$PROJECT_BASE/01.RawData"
QUALITYDIR_F="$PROJECT_BASE/fastqc_output/forward_reads"
QUALITYDIR_R="$PROJECT_BASE/fastqc_output/reverse_reads"

cd "$PROJECT_BASE"

mkdir -p "$DATADIR" "$QUALITYDIR_F" "$QUALITYDIR_R" slurm_logs
cd "$DATADIR"

# Download MicrobesNG url list (edit to your project URL)
wget -nc https://microbesng-data.s3-eu-west-1.amazonaws.com/projects/.../untrimmed_urls.txt

# Download each URL listed
if [ -s untrimmed_urls.txt ]; then
  while read -r url; do
    echo "Downloading: $url"
    wget -nc "$url"
  done < untrimmed_urls.txt
else
  echo "Warning: untrimmed_urls.txt missing or empty"
fi

# Create names.txt (overwrite each run)
: > names.txt # Empty names.txt 
for file in *_1.fastq.gz; do
  [ -e "$file" ] || continue
  base=$(basename "$file" _1.fastq.gz)
  echo "$base" >> names.txt
done

# Run FastQC
while read -r sample; do
  echo "Running FastQC on forward read: ${sample}_1.fastq.gz"
  fastqc "${sample}_1.fastq.gz" -o "$QUALITYDIR_F"
done < names.txt

while read -r sample; do
  echo "Running FastQC on reverse read: ${sample}_2.fastq.gz"
  fastqc "${sample}_2.fastq.gz" -o "$QUALITYDIR_R"
done < names.txt

# MultiQC summary
cd "$QUALITYDIR_F"
multiqc . -o . || true
cd "$QUALITYDIR_R"
multiqc . -o . || true

echo "FastQC and MultiQC analyses complete. Reports in: $QUALITYDIR_F and $QUALITYDIR_R"
