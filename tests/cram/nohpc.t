  $ IN=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/index_default.mmi --log-level=DEBUG 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 1 (glob)
  $ $__PBTEST_PBMM2_EXE align $IN $CRAMTMP/index_default.mmi $CRAMTMP/out_hpc.bam --log-level=DEBUG 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 1 (glob)
  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/index_default.mmi --log-level=DEBUG --no-hpc 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 0 (glob)
  $ $__PBTEST_PBMM2_EXE align $IN $CRAMTMP/index_default.mmi $CRAMTMP/out_hpc.bam --log-level=DEBUG --no-hpc 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 0 (glob)
  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/index_default.mmi --log-level=DEBUG --preset ISOSEQ 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 0 (glob)
  $ $__PBTEST_PBMM2_EXE align $IN $CRAMTMP/index_default.mmi $CRAMTMP/out_hpc.bam --log-level=DEBUG --preset ISOSEQ 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 0 (glob)
  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/index_default.mmi --log-level=DEBUG --preset CCS 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 0 (glob)
  $ $__PBTEST_PBMM2_EXE align $IN $CRAMTMP/index_default.mmi $CRAMTMP/out_hpc.bam --log-level=DEBUG --preset CCS 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 0 (glob)

  $ samtools view -H $TESTDIR/data/median.bam > $CRAMTMP/single_subread.sam
  $ samtools view -h $TESTDIR/data/median.bam | head -n 1 >> $CRAMTMP/single_subread.sam
  $ samtools view -bS $CRAMTMP/single_subread.sam > $CRAMTMP/single_subread.bam
  $ IN=$CRAMTMP/single_subread.bam

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/out_hpc.bam --log-level=DEBUG 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 1 (glob)
  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/out_hpc.bam --log-level=DEBUG --no-hpc 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 0 (glob)
  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/out_hpc.bam --log-level=DEBUG --preset ISOSEQ 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 0 (glob)
  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/out_hpc.bam --log-level=DEBUG --preset CCS 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : 0 (glob)
