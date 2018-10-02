  $ IN=$TESTDIR/data/m54019_171011_032401_tiny.subreadset.xml
  $ REF=$TESTDIR/data/lambdaNEB_BsaAI_allFrags_incLeftRightEnds_unrolled_250k.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/zmw.alignmentset.xml --log-level INFO --zmw 2>&1| grep INFO | grep "mapped read length"
  *Max mapped read length : 156345 (glob)
  *Mean mapped read length : 125536 (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/hqregion.alignmentset.xml --log-level INFO --hqregion 2>&1| grep INFO | grep "mapped read length"
  *Max mapped read length : 125134 (glob)
  *Mean mapped read length : 82758 (glob)
