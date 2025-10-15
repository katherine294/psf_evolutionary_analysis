Run scripts 1 - 8.

# itol was used to create the final tree figure.

SBATCH 01.download_fastqc.sh

SBATCH 02.trimmomatic.sh

SBATCH 03.variant_calling.sh

SBATCH 04.core_alignment.sh

SBATCH 05.recombination_detection.sh

SBATCH 06.MGE_masking.sh

SBATCH 07.merge_sort_masking_bed.sh

SBATCH 08.Masked_core_SNP_tree.sh
