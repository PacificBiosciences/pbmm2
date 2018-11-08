  $ BAM=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ samtools view ${BAM} | awk '{ print ">"$1"\n"$10 }' > $CRAMTMP/median.fasta
  $ FASTA=$CRAMTMP/median.fasta
  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/ecoli.mmi --log-level FATAL
  $ REF=$CRAMTMP/ecoli.mmi 

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/rg.bam --log-level FATAL
  $ samtools view -H $CRAMTMP/rg.bam | grep "^@RG"
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)
  $ rm $CRAMTMP/rg.bam

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/rg.bam --log-level FATAL --rg '@RG'
  *Invalid @RG line. Missing ID field. Please provide following format: '@RG\tID:xyz\tSM:abc'* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/rg.bam --log-level FATAL --rg '@XY'
  *Invalid @RG line. Missing ID field. Please provide following format: '@RG\tID:xyz\tSM:abc'* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/rg.bam --log-level FATAL --sample newSampleName
  $ samtools view -H $CRAMTMP/rg.bam | grep "^@RG"
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:newSampleName\tPM:SEQUEL (esc)
  $ rm $CRAMTMP/rg.bam

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/rg.bam --log-level FATAL --rg '@RG\tID:myid'
  $ samtools view -H $CRAMTMP/rg.bam | grep "^@RG"
  @RG\tID:myid\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)
  $ rm $CRAMTMP/rg.bam

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/rg.bam --log-level FATAL --rg '@RG\tID:myid\tSM:mysample'
  $ samtools view -H $CRAMTMP/rg.bam | grep "^@RG"
  @RG\tID:myid\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:mysample\tPM:SEQUEL (esc)
  $ rm $CRAMTMP/rg.bam

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/rg.bam --log-level FATAL --rg '@RG\tID:myid\tSM:mysample' --sample newSampleName
  $ samtools view -H $CRAMTMP/rg.bam | grep "^@RG"
  @RG\tID:myid\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:newSampleName\tPM:SEQUEL (esc)
  $ rm $CRAMTMP/rg.bam

  $ $__PBTEST_PBMM2_EXE align $REF $BAM $CRAMTMP/rg.bam --log-level FATAL --rg '@RG\tID:myid'
  *Cannot override read groups with BAM input. Remove option --rg.* (glob)
  [1]
