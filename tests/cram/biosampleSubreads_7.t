  $ MERGED="$TESTDIR"/data/merged.dataset.xml
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align "$MERGED" "$REF" "$CRAMTMP"/out7.bam --sample testSample
  $ "$SAMTOOLS" view -H "$CRAMTMP"/out7.bam | grep "@RG"
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
