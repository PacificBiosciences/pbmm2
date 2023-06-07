  $ IN="$TESTDIR"/data/inv.fasta
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/long.bam --preset SUBREAD --sort --log-level FATAL
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/short.bam --preset SUBREAD --sort --short-sa-cigar --log-level FATAL

  $ samtools view "$CRAMTMP"/long.bam | perl -ne '/.*SA:Z:([^\t]*).*/ && print($1)' | tr '=' "M" > "$CRAMTMP"/long.sa
  $ samtools view "$CRAMTMP"/short.bam | perl -ne '/.*SA:Z:([^\t]*).*/ && print($1)' > "$CRAMTMP"/short.sa
  $ diff "$CRAMTMP"/short.sa "$CRAMTMP"/long.sa
