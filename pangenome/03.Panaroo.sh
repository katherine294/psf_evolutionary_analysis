#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=6750M
#SBATCH --ntasks=54
#SBATCH --nodes=1
#SBATCH --output=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/panaroo_%A_%a.out
#SBATCH --error=/rds/homes/k/kgh742/psf_wgs_project/slurm_logs/panaroo_%A_%a.err

module purge
module load bear-apps/2022a
module load panaroo/1.5.0-foss-2022a



mkdir pangenome_gff_input

cp */*.gff3 */*.fna pangenome_gff_input

# Create panaroo input list

ls *.gff3 | while read gff; do
    base=$(basename "$gff" .gff3)
    echo "pangenome_gff_input/${base}.gff3  pangenome_gff_input/${base}.fna"
done > psf_reduced_500.txt


# Run Panaroo
panaroo -i panaroo_input-psav.txt -o panaroo_merged_psav-1704 \
-t 50 \
--clean-mode strict \
-c 0.98 \
-f 0.7 \
--len_dif_percent 0.98 \
--aligner mafft \
--core_threshold 0.98 \
-a pan \
--merge_paralogs
