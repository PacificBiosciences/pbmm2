  $ IN=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/aligned.bam --cstag
  $ samtools view $CRAMTMP/aligned.bam | head -n 1 | grep -c '\tcs:'
  1