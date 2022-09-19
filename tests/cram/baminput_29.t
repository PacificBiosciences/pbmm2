  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/bestn1.bam --best-n 1
  $ samtools view "$CRAMTMP"/bestn1.bam | wc -l | tr -d ' '
  52
