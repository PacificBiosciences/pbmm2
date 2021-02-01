  $ IN="$TESTDIR"/data/median.bam
  $ MERGED="$TESTDIR"/data/mergd.dataset.xml
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align "$IN" "$REF" "$CRAMTMP"/out2.bam --sample testSample
  $ "$SAMTOOLS" view -H "$CRAMTMP"/out2.bam | grep "@RG"
  *\tSM:testSample\t* (glob)
