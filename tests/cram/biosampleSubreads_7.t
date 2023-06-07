  $ MERGED="$TESTDIR"/data/merged.dataset.xml
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$MERGED" "$REF" "$CRAMTMP"/out7.bam --preset SUBREAD --sample testSample
  $ samtools view -H "$CRAMTMP"/out7.bam | grep "@RG"
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
