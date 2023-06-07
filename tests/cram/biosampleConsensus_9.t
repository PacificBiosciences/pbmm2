  $ MERGED="$TESTDIR"/data/merged.consensusreadset.xml
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$MERGED" "$REF" "$CRAMTMP"/biosample_consensus_9.bam --preset CCS
  *Offending bio sample names. BAM contains 'bamSample' and XML contains 'UCLA 1023'. Will ignore XML bio sample name.* (glob)
  $ samtools view -H "$CRAMTMP"/biosample_consensus_9.bam | grep "@RG"
  *\tSM:bamSample\t* (glob)
  *\tSM:test test\t* (glob)
  *\tSM:3260208_188nM-GTAC_4xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494\t* (glob)
