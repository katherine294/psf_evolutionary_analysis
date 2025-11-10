Run scripts 1 - 8.

# itol was used to create the final tree figure.

sbatch 01.download_fastqc.sh

sbatch 02.trimmomatic.sh

sbatch 03.variant_calling.sh

sbatch 04.core_alignment.sh

sbatch 05.recombination_detection.sh

sbatch 06.MGE_masking.sh

sbatch 07.merge_sort_masking_bed.sh

sbatch 08.Masked_core_SNP_tree.sh
