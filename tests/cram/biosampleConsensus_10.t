  $ MERGED="$TESTDIR"/data/merged.consensusreadset.xml
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$MERGED" "$REF" "$CRAMTMP"/biosample_consensus_10.bam --preset CCS --sample testSample
  *Offending bio sample names. BAM contains 'bamSample' and XML contains 'UCLA 1023'. Will ignore XML bio sample name.* (glob)
  $ samtools view -H "$CRAMTMP"/biosample_consensus_10.bam | grep "@RG"
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)

