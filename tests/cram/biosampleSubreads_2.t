  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/out2.bam --preset SUBREAD --sample testSample
  $ samtools view -H "$CRAMTMP"/out2.bam | grep "@RG"
  *\tSM:testSample\t* (glob)
