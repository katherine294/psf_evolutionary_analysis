#!/usr/bin/env python3
import pandas as pd
import glob
import os

# Path to your summary files (current directory by default)
summary_files = glob.glob("*.mosdepth.summary.txt")

if not summary_files:
    raise FileNotFoundError("No *.mosdepth.summary.txt files found in this directory")

data = []

for file in summary_files:
    strain = os.path.basename(file).split(".")[0]  # Get strain name from filename
    df = pd.read_csv(file, sep="\t")

    # Filter only main contig lines (exclude "_region" and "total")
    df_filtered = df[~df['chrom'].str.contains("_region|total")]

    # Extract mean coverages
    mean_cov = df_filtered.set_index("chrom")["mean"].to_dict()

    # Add strain name at start
    row = {"Strain": strain}
    row.update(mean_cov)
    data.append(row)

# Create DataFrame
coverage_df = pd.DataFrame(data)

# Optional: reorder columns
coverage_df = coverage_df[["Strain"] + sorted([col for col in coverage_df.columns if col != "Strain"])]

# Save to file with date
output_file = "coverage_summary_{}.tsv".format(pd.Timestamp.now().strftime("%y%m%d"))
coverage_df.to_csv(output_file, sep="\t", index=False)

print(f"Coverage summary written to: {output_file}")
print(coverage_df.head())
