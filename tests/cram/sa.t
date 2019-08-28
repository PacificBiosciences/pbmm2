  $ IN=$TESTDIR/data/bnd.bam
  $ IN_FASTA=$TESTDIR/data/bnd.fasta
  $ REF=$TESTDIR/data/bnd-ref.fasta

  $ minimap2 -x map-pb -a --eqx -L -O 5,56 -E 4,1 -B 5 --secondary=no -z 400,50 -r 2k -Y ${REF} ${IN_FASTA} 2>/dev/null | samtools sort > $CRAMTMP/mm2.bam
  $ $__PBTEST_PBMM2_EXE align ${IN} ${REF} $CRAMTMP/pbmm2.bam --sort -j 1 -J 1
  $ samtools view $CRAMTMP/mm2.bam | rev | cut -f 1 | rev > mm2.sa
  $ samtools view $CRAMTMP/pbmm2.bam | rev | cut -f 1 | rev > pbmm2.sa
  $ diff mm2.sa pbmm2.sa
