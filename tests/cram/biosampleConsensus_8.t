  $ IN2="$TESTDIR"/data/m54075_180905_225130.ccs.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN2" "$REF" "$CRAMTMP"/biosample_consensus_8.bam --sample ""
  $ "$SAMTOOLS" view -H "$CRAMTMP"/biosample_consensus_8.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)
