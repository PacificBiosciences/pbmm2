  $ IN=$TESTDIR/data/median.bam
  $ MERGED=$TESTDIR/data/merged.dataset.xml
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/out.bam
  $ samtools view -H $CRAMTMP/out.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/out2.bam --sample-name testSample
  $ samtools view -H $CRAMTMP/out2.bam | grep "@RG"
  *\tSM:testSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/out3.bam --sample-name "   TEST bla   "
  $ samtools view -H $CRAMTMP/out3.bam | grep "@RG"
  *\tSM:TEST_bla\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/out4.bam --sample-name "      "
  $ samtools view -H $CRAMTMP/out4.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/out5.bam --sample-name ""
  $ samtools view -H $CRAMTMP/out5.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $MERGED $REF $CRAMTMP/out6.bam
  $ samtools view -H $CRAMTMP/out6.bam | grep "@RG"
  *\tSM:3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494\t* (glob)
  *\tSM:test_test\t* (glob)
  *\tSM:UCLA_1023\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $MERGED $REF $CRAMTMP/out7.bam --sample-name testSample
  $ samtools view -H $CRAMTMP/out7.bam | grep "@RG"
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
