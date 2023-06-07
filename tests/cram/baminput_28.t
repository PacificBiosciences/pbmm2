  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

Test that median filter does not fail
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/median_output.bam --preset SUBREAD --median-filter
