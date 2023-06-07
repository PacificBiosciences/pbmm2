  $ IN2="$TESTDIR"/data/m54075_180905_225130.ccs.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN2" "$REF" "$CRAMTMP"/biosample_consensus_2.bam --preset CCS
  $ samtools view -H "$CRAMTMP"/biosample_consensus_2.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)
