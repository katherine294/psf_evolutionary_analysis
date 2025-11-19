#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <blast_outfmt6.txt> <output.bed>"
  exit 1
fi

INPUT="$1"
OUTPUT="$2"

echo "# chrom  start  end  .  .  strand" > "$OUTPUT"

awk '
# Skip comment/blank lines
/^[[:space:]]*#/ || NF==0 { next }

{
    sseqid = $2      # chrom
    sstart = $9      # subject start
    send   = $10     # subject end

    # Determine coordinates and strand
    if (sstart <= send) {
        start = sstart - 1
        end   = send
        strand = "+"
    } else {
        start = send - 1
        end   = sstart
        strand = "-"
    }

    # Print chrom, start, end, empty name, empty score, strand
    print sseqid, start, end, ".", ".", strand
}
' OFS="\t" "$INPUT" >> "$OUTPUT"

echo "Done. BED written to $OUTPUT"
