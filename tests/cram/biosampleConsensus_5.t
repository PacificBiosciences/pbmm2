  $ IN="$TESTDIR"/data/m54075_180905_221350.ccs.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align "$IN" "$REF" "$CRAMTMP"/biosample_consensus_5.bam --sample "      "
  $ "$SAMTOOLS" view -H "$CRAMTMP"/biosample_consensus_5.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)
