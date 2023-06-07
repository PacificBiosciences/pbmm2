  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta
  $ source "$TESTDIR"/sorttags.sh

  $ samtools view -h "$IN" | head -n 20 | samtools view -bS > "$CRAMTMP"/median_12.bam
  $ IN="$CRAMTMP"/median_12.bam
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/sorted.bam --preset SUBREAD --sort --log-level FATAL
  $ "$PBMM2" align -j 1 "$CRAMTMP"/sorted.bam "$REF" "$CRAMTMP"/idempotent.bam --preset SUBREAD --sort --log-level FATAL
  $ sorttags "$CRAMTMP"/sorted.bam | sort > "$CRAMTMP"/sorted.sam
  $ sorttags "$CRAMTMP"/idempotent.bam | sort > "$CRAMTMP"/idempotent.sam
  $ diff "$CRAMTMP"/sorted.sam "$CRAMTMP"/idempotent.sam
