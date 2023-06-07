  $ IN="$TESTDIR"/data/m54075_180905_221350.ccs.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/biosample_consensus_3.bam --preset CCS --sample testSample
  $ samtools view -H "$CRAMTMP"/biosample_consensus_3.bam | grep "@RG"
  *\tSM:testSample\t* (glob)
