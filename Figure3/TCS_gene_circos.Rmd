---
title: "tcs_circos"
author: "katherine hinton"
date: "2026-02-20"
output: html_document
---

----------------- Filter strain vcf files to retain mutations in specific genes ---------------------
 
1. Convert vcf files output from snippy to bedfile
```{bash}
#!/bin/bash
module load bear-apps/2023a
module load BCFtools/1.20-GCC-12.3.0


# Define the directory where VCF files will be copied
output_dir="circos_strains_vcf"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"


while read -r file; do
    if [[ -f "${file}.vcf" ]]; then
        # Extract SNP positions, remove header, and ensure format is correct for Circos
        bcftools query -f '%CHROM\t%POS\t%POS\t[%INFO/TYPE]\n' "${file}.vcf" | grep -v "^#" > "$output_dir/${file}_circos.txt"
        echo "Processed: ${file}.vcf -> $output_dir/${file}_circos.txt"
    else
        echo "Warning: ${file}.vcf not found, skipping..."
    fi
done < strains.txt
```

2. Prepare txt list of gene IDs and filter reference gff with these genes
```{bash}
#!/bin/bash

# Files
gff="246539E_W163A3B1.gff"
gene_list="TCS_gene_list.txt" # list of TCS genes output from KEGG enrcihment analysis
bed_output="TCS_genes.bed" # Only keep TCS genes as bed file 

# Make sure output file is empty first
> TCS_genes.bed

# Loop over each gene in gene.txt
while read -r GENE; do
        GENE=$(echo "$GENE" | tr -d '\r')   # remove carriage return
        awk -v gene="$GENE" '$0 ~ gene {print $1, $4-1, $5, gene}' "$gff" >> "$bed_output"
done < "$gene_list"
```

3. Filter strain bedfiles to only keep TCS mutations using filtered gff reference file
```{bash}
#!/bin/bash

module purge
module load bear-apps/2024a
module load BEDTools/2.31.1-GCC-13.3.0

INPATH="/rds/homes/k/kgh742/psf_wgs_project/01.RawData/02.TrimmedReads/03.Snippy_W163_chrom/vcffiles/circos_strains_vcf"
OUTPATH="/rds/homes/k/kgh742/psf_wgs_project/01.RawData/02.TrimmedReads/03.Snippy_W163_chrom/vcffiles/circos_strains_vcf_filtered"
STRAIN_FILE="/rds/homes/k/kgh742/psf_wgs_project/01.RawData/02.TrimmedReads/03.Snippy_W163_chrom/vcffiles/strains.txt"
TCS_BEDFILE="/rds/homes/k/kgh742/psf_wgs_project/01.RawData/02.TrimmedReads/03.Snippy_W163_chrom/vcffiles/TCS_genes_tab.bed"

# Make sure output directory exists
mkdir -p "$OUTPATH"

while read -r STRAIN; do
    INFILE="${INPATH}/${STRAIN}.bed"
    OUTFILE="${OUTPATH}/${STRAIN}_filtered_variants.bed"

    # Intersect
    bedtools intersect -a "$INFILE" -b "$TCS_BEDFILE" -wa > "$OUTFILE"
done < "$STRAIN_FILE"

```

4. Convert mutation data to circos format including colour metadata 
```{bash}
#!/bin/bash

# Ensure your input text file with strain names is 'strains.txt'

input_dir="circos_strains_vcf_filtered"
output_dir="circos_modified_filtered_tcs"

mkdir -p "$output_dir"

# Loop through each strain name in strains.txt
while read -r strain; do
    echo "Processing strain: $strain"  # Debugging statement to check reading of strains.txt

# Construct the filename
    file="${input_dir}/${strain}_filtered_variants.bed"

    # Check if the file exists
    if [[ -f "$file" ]]; then
        # Create a new output file to store the modified results
        output="$output_dir/${strain}_modified_circos.txt"

        # Process each line and replace mutation type with corresponding fill_color
        while read -r line; do
            if [[ $line =~ "snp" ]]; then
                echo "$line" | sed 's/snp/fill_color=vdred,stroke_color=vdred,stroke_thickness=2/' >> "$output"
            elif [[ $line =~ "complex" ]]; then
                echo "$line" | sed 's/complex/fill_color=vvdblue,stroke_color=vvdblue,stroke_thickness=2/' >> "$output"
            elif [[ $line =~ "ins" ]]; then
                echo "$line" | sed 's/ins/fill_color=vdpurple,stroke_color=vdpurple,stroke_thickness=2/' >> "$output"
            elif [[ $line =~ "del" ]]; then
                echo "$line" | sed 's/del/fill_color=vdgreen,stroke_color=vdgreen,stroke_thickness=2/' >> "$output"
            elif [[ $line =~ "mnp" ]]; then
                echo "$line" | sed 's/mnp/fill_color=orange,stroke_color=orange,stroke_thickness=2/' >> "$output"

            else
                # If no match, just echo the line as is
                echo "$line" >> "$output"
            fi
        done < "$file"
```

5. Append the radius position information for each strain so they fit onto the circos plot from the outside into the centre
```{bash}
# 1. Look at your coresnps.txt file and order strain according to most to least number of snps
# 2. Save strain order in new txt file called strains_ordered.txt
# 3. Choose start and end point on circos plot - i am using 0.92r to 0r. I have 123 strains. So, each strain has 0.92/123 of space (~0.075)

#!/bin/bash

# Set starting r0 and step size
r0=0.93         # This is your outermost point 
step=0.0075     # Calculate this value by r0/strain count

while read strain; do
    r1=$(echo "$r0 - $step" | bc)

    # Format r0 and r1 with trailing "r"
    r0_fmt=$(printf "%.5f" "$r0")r
    r1_fmt=$(printf "%.5f" "$r1")r

    bedfile="${strain}_modified_circos.txt"
    outfile="${strain}_updated_all.txt"

    if [[ -f "$bedfile" ]]; then
        awk -v r0="$r0_fmt" -v r1="$r1_fmt" '{
            if ($0 ~ /^#/ || NF==0) { print; next }
            $4 = $4 ",r0=" r0 ",r1=" r1
            print
        }' "$bedfile" > "$outfile"
    else
        echo "WARNING: $bedfile not found, skipping"
    fi

    # Move r0 for next strain
    r0=$r1

done < strains_all_ordered.txt

# Merge step
cat *_updated_all.txt > merged_all_circos.txt

echo "Finished! Output: merged_for_circos.txt"

```

6. Prepare reference for circos by seperating forward and reverse CDS and converting reference to kareotype file 

Extract CDS info from reference gff file
```{bash}

 awk '$3 == "CDS" {print $1"\t"$4"\t"$5"\t"$9"\t"$7}' W163a3b1_genes.gff3 > cds.bed
```

Extract forward and reverse reads into individual bed files
```{python}
# Open the input CDS file
input_file = "cds.bed"
forward_file = "cds_forward.bed"
reverse_file = "cds_reverse.bed"

# Open the output files for the forward and reverse CDS
with open(forward_file, 'w') as fwd, open(reverse_file, 'w') as rev:
    with open(input_file, 'r') as infile:
        for line in infile:
            columns = line.strip().split('\t')
            # Extract coordinates (start, end) from the CDS entry
            chromosome = columns[0]
            start = columns[1]
            end = columns[2]
            
            # Check the strand based on the direction (+ or -)
            if columns[4] == '+':
                # For forward strand, write to the forward file
                fwd.write(f"{chromosome}\t{start}\t{end}\n")
            elif columns[4] == '-':
                # For reverse strand, write to the reverse file
                rev.write(f"{chromosome}\t{start}\t{end}\n")

```

Make .fai file from fasta
```{bash}
module load bear-apps/2023a
module load SAMtools/1.21-GCC-12.3.0

samtools faidx 246539E_W163A3B1.fasta

```

Make kareotype file from .fai file
```{bash}
awk '{print "chr - "$1" "$1" 0 "$2" black"}' 246539E_W163A3B1.fasta.fai > chr.kar
```

-------------- Circo config file ----------------------------

```{circos}

karyotype = chr.kar
chromosomes_units = 1000000

<<include colors_fonts_patterns.conf>>
<<include housekeeping.conf>>

# IMAGE
<image>
<<include image.conf>>
</image>

# IDEOGRAM
<ideogram>
<spacing>
default = 0.0005r
</spacing>
radius           = 0.8r
thickness        = 5p
fill             = yes
stroke_color     = black
stroke_thickness = 5p
#show_label       = yes
#label_radius     = 1.1r
#label_size       = 80p
#label_parallel   = yes
show_background = yes

<backgrounds>
<background>
color = vvlgrey
# the "r" suffix indicates position relative to track data range
y0    = 0.93r
y1    = 0.0r
</background>
</backgrounds>

</ideogram>

<highlights>
z = 0
fill_color = green

<highlight>
file       = cds_forward_1.bed
r0         = 0.95r
r1         = 1.0r
fill_color = vvdorange
stroke_color = vvdorange
stroke_thickness = 2
</highlight>

<highlight>
file       = cds_reverse_1.bed
r0         = 1.0r
r1         = 1.05r
fill_color = ldorange
stroke_color = ldorange
stroke_thickness = 2
</highlight>

<highlight>
file = tcs_merged_circos_v2.txt
</highlight>

</highlights>


<<include ticks.conf>>

```  
