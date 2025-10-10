#### GENERATE AN UNMASKED CORE ALIGNMENT USING SNIPPY-CORE, FOR RECOMBINATION DETECTION ####

ALN_PREFIX="W163a3b1_core_genome_unmasked"
STRAINS=$(cat Names.txt)

cd $SNIPPY_OUTPUT

snippy-core \
  --prefix "$ALN_PREFIX" \
  --ref "$REF_DIR/246539E_W163A3B1_chrom.fasta" \
  "$STRAINS"

ALN_CAND="$SNIPPY_OUTPUT/${ALN_PREFIX}.full.aln"

# Clean alignment (replace ambiguous or low-quality sites with "X")
snippy-clean_full_aln $ALN_CAND > $SNIPPY_OUTPUT/${ALN_PREFIX}_clean.full.aln
echo "Cleaned alignment written to: ${SNIPPY_OUTPUT}/${ALN_PREFIX}_clean.full.aln"

echo "Per-sample snippy complete"
