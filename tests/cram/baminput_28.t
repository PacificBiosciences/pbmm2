  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

Test that median filter does not fail
  $ "$PBMM2" align "$IN" "$REF" "$CRAMTMP"/median_output.bam --median-filter
