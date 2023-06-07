  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/out.bam --preset SUBREAD
  $ samtools view -H "$CRAMTMP"/out.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)
