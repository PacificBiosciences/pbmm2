  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/out3.bam --preset SUBREAD --sample "   TEST bla   "
  $ samtools view -H "$CRAMTMP"/out3.bam | grep "@RG"
  *\tSM:TEST bla\t* (glob)
