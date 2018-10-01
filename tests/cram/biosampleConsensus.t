  $ IN=$TESTDIR/data/m54075_180905_221350.ccs.bam
  $ MERGED=$TESTDIR/data/merged.consensusreadset.xml
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs.bam
  $ samtools view -H $CRAMTMP/ccs.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs2.bam --sample-name testSample
  $ samtools view -H $CRAMTMP/ccs2.bam | grep "@RG"
  *\tSM:testSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs3.bam --sample-name "   TEST bla   "
  $ samtools view -H $CRAMTMP/ccs3.bam | grep "@RG"
  *\tSM:TEST_bla\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs4.bam --sample-name "      "
  $ samtools view -H $CRAMTMP/ccs4.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs5.bam --sample-name ""
  $ samtools view -H $CRAMTMP/ccs5.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $MERGED $REF $CRAMTMP/ccs6.bam
  $ samtools view -H $CRAMTMP/ccs6.bam | grep "@RG"
  *\tSM:UCLA_1023\t* (glob)
  *\tSM:test_test\t* (glob)
  *\tSM:3260208_188nM-GTAC_4xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $MERGED $REF $CRAMTMP/ccs7.bam --sample-name testSample
  $ samtools view -H $CRAMTMP/ccs7.bam | grep "@RG"
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
