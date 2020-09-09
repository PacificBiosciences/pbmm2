  $ IN=$TESTDIR/data/m54075_180905_221350.ccs.bam
  $ IN2=$TESTDIR/data/m54075_180905_225130.ccs.bam
  $ MERGED=$TESTDIR/data/merged.consensusreadset.xml
  $ NO_SM_BIOSAMPLES=$TESTDIR/data/no_sm_biosamples.consensusreadset.xml
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs.bam
  $ samtools view -H $CRAMTMP/ccs.bam | grep "@RG"
  *\tSM:bamSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN2 $REF $CRAMTMP/ccs_2.bam
  $ samtools view -H $CRAMTMP/ccs_2.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs2.bam --sample testSample
  $ samtools view -H $CRAMTMP/ccs2.bam | grep "@RG"
  *\tSM:testSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs3.bam --sample "   TEST bla   "
  $ samtools view -H $CRAMTMP/ccs3.bam | grep "@RG"
  *\tSM:TEST_bla\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs4.bam --sample "      "
  $ samtools view -H $CRAMTMP/ccs4.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ccs5.bam --sample ""
  $ samtools view -H $CRAMTMP/ccs5.bam | grep "@RG"
  *\tSM:bamSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN2 $REF $CRAMTMP/ccs4_2.bam --sample "      "
  $ samtools view -H $CRAMTMP/ccs4_2.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN2 $REF $CRAMTMP/ccs5_2.bam --sample ""
  $ samtools view -H $CRAMTMP/ccs5_2.bam | grep "@RG"
  *\tSM:UnnamedSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $MERGED $REF $CRAMTMP/ccs6.bam
  $ samtools view -H $CRAMTMP/ccs6.bam | grep "@RG"
  *\tSM:bamSample\t* (glob)
  *\tSM:test_test\t* (glob)
  *\tSM:test_test\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $MERGED $REF $CRAMTMP/ccs7.bam --sample testSample
  $ samtools view -H $CRAMTMP/ccs7.bam | grep "@RG"
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)
  *\tSM:testSample\t* (glob)

  $ $__PBTEST_PBMM2_EXE align $NO_SM_BIOSAMPLES $REF $CRAMTMP/ccs8.bam
  *<BioSamples> list element is present in dataset XML, but SM tags are missing* (glob)
  [1]
