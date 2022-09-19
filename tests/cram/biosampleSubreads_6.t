  $ MERGED="$TESTDIR"/data/merged.dataset.xml
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$MERGED" "$REF" "$CRAMTMP"/out6.bam
  $ samtools view -H "$CRAMTMP"/out6.bam | grep "@RG"
  *\tSM:3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494\t* (glob)
  *\tSM:test_test\t* (glob)
  *\tSM:UCLA_1023\t* (glob)
