  $ IN="$TESTDIR"/data/m54075_180905_221350.ccs.bam
  $ MERGED="$TESTDIR"/data/merged.consensusreadset.xml
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/biosample_consensus_1.bam
  $ samtools view -H "$CRAMTMP"/biosample_consensus_1.bam | grep "@RG"
  *\tSM:bamSample\t* (glob)
