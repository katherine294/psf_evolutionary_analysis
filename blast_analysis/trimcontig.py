#!/usr/bin/env python

from Bio import SeqIO
import sys

# Input and output files as command-line arguments
input_file = sys.argv[1]  # The input FASTA file
output_file = sys.argv[2]  # The output filtered FASTA file

# Minimum contig length
min_length = 500

# Filter sequences longer than 500 bp
with open(output_file, "w") as out_fasta:
    for record in SeqIO.parse(input_file, "fasta"):
        if len(record.seq) >= min_length:
            SeqIO.write(record, out_fasta, "fasta")
