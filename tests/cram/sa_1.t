  $ IN="$TESTDIR"/data/bnd.bam
  $ IN_FASTA="$TESTDIR"/data/bnd.fasta
  $ REF="$TESTDIR"/data/bnd-ref.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/pbmm2.bam --preset SUBREAD --sort -j 1 -J 1 --short-sa-cigar

  $ samtools view "$CRAMTMP"/pbmm2.bam | sed -E 's/^.*(SA:Z:[ !-~]*).*/\1/g'
  SA:Z:bnd_ref_V,601,+,600S600M,60,0;
  SA:Z:bnd_ref_U,601,+,600M600S,60,0;
  SA:Z:bnd_ref_Y,601,-,600M600S,60,0;
  SA:Z:bnd_ref_W,601,+,600M600S,60,0;
  SA:Z:bnd_ref_Z,601,+,600S600M,60,0;
  SA:Z:bnd_ref_X,601,-,600S600M,60,0;
