  $ IN=$TESTDIR/data/bnd.bam
  $ IN_FASTA=$TESTDIR/data/bnd.fasta
  $ REF=$TESTDIR/data/bnd-ref.fasta

  $ minimap2 -x map-pb -a --eqx -L -O 5,56 -E 4,1 -B 5 --secondary=no -z 400,50 -r 2k -Y ${REF} ${IN_FASTA} 2>/dev/null | samtools sort > $CRAMTMP/mm2.bam
  $ $__PBTEST_PBMM2_EXE align ${IN} ${REF} $CRAMTMP/pbmm2.bam --sort -j 1 -J 1
  $ samtools view $CRAMTMP/mm2.bam | rev | cut -f 1 | rev > mm2.sa
  $ samtools view $CRAMTMP/pbmm2.bam | rev | cut -f 1 | rev > pbmm2.sa
  $ diff mm2.sa pbmm2.sa > $CRAMTMP/idem.diff
  [1]
  $ cat $CRAMTMP/idem.diff | tr -d '>' | tr -d '<'
  1,6c1,6
   SA:Z:bnd_ref_V,601,+,600S600M,60,0;
   SA:Z:bnd_ref_U,601,+,600M600S,60,0;
   SA:Z:bnd_ref_Y,601,-,600M600S,60,0;
   SA:Z:bnd_ref_W,601,+,600M600S,60,0;
   SA:Z:bnd_ref_Z,601,+,600S600M,60,0;
   SA:Z:bnd_ref_X,601,-,600S600M,60,0;
  ---
   SA:Z:bnd_ref_V,601,+,600S599M1S,60,0;
   SA:Z:bnd_ref_U,601,+,599M601S,60,0;
   SA:Z:bnd_ref_Y,602,-,1S599M600S,60,0;
   SA:Z:bnd_ref_W,601,+,599M601S,60,0;
   SA:Z:bnd_ref_Z,601,+,600S599M1S,60,0;
   SA:Z:bnd_ref_X,602,-,601S599M,60,0;
