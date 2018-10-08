  $ IN=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted.bam --sort --log-level FATAL
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/sorted.bam $REF $CRAMTMP/idempotent.bam --sort --log-level FATAL
  $ samtools view $CRAMTMP/sorted.bam | sort -k1,1g -k3,3g -k4,4g > $CRAMTMP/sorted.sam
  $ samtools view $CRAMTMP/idempotent.bam | sort -k1,1g -k3,3g -k4,4g > $CRAMTMP/idempotent.sam
  $ diff $CRAMTMP/sorted.sam $CRAMTMP/idempotent.sam

  $ rm $CRAMTMP/sorted.sam $CRAMTMP/idempotent.sam $CRAMTMP/sorted.bam $CRAMTMP/idempotent.bam

  $ IN=$TESTDIR/data/m54075_180905_221350.ccs.bam
  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted.bam --sort --log-level FATAL
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/sorted.bam $REF $CRAMTMP/idempotent.bam --sort --log-level FATAL
  $ samtools view $CRAMTMP/sorted.bam | sort -k1,1g -k3,3g -k4,4g > $CRAMTMP/sorted.sam
  $ samtools view $CRAMTMP/idempotent.bam | sort -k1,1g -k3,3g -k4,4g > $CRAMTMP/idempotent.sam
  $ diff $CRAMTMP/sorted.sam $CRAMTMP/idempotent.sam

  $ rm $CRAMTMP/sorted.sam $CRAMTMP/idempotent.sam $CRAMTMP/sorted.bam $CRAMTMP/idempotent.bam

  $ IN=$TESTDIR/data/m54019_171011_032401_tiny.subreadset.xml
  $ REF=$TESTDIR/data/lambdaNEB_BsaAI_allFrags_incLeftRightEnds_unrolled_250k.fasta
  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted.bam --sort --log-level FATAL --zmw
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/sorted.bam $REF $CRAMTMP/idempotent.bam --sort --zmw --log-level FATAL
  $ samtools view $CRAMTMP/sorted.bam | sort -k1,1g -k3,3g -k4,4g > $CRAMTMP/sorted.sam
  $ samtools view $CRAMTMP/idempotent.bam | sort -k1,1g -k3,3g -k4,4g > $CRAMTMP/idempotent.sam
  $ diff $CRAMTMP/sorted.sam $CRAMTMP/idempotent.sam
