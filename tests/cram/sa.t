  $ IN="$TESTDIR"/data/bnd.bam
  $ IN_FASTA="$TESTDIR"/data/bnd.fasta
  $ REF="$TESTDIR"/data/bnd-ref.fasta

  $ "$MINIMAP2" -x map-pb -a --eqx -L -O 5,56 -E 4,1 -B 5 --secondary=no -z 400,50 -r 2k -Y "$REF" "$IN_FASTA" 2>/dev/null | "$SAMTOOLS" sort > "$CRAMTMP"/mm2.bam
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/pbmm2.bam --sort -j 1 -J 1
  $ "$SAMTOOLS" view "$CRAMTMP"/mm2.bam | rev | cut -f 1 | rev > mm2.sa
  $ "$SAMTOOLS" view "$CRAMTMP"/pbmm2.bam | rev | cut -f 1 | rev > pbmm2.sa
  $ diff mm2.sa pbmm2.sa
