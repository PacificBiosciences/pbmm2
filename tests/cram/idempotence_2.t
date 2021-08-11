  $ IN="$TESTDIR"/data/m54075_180905_221350.ccs.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta
  $ source "$TESTDIR"/sorttags.sh

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/sorted.bam --sort --log-level FATAL
  $ "$PBMM2" align -j 1 "$CRAMTMP"/sorted.bam "$REF" "$CRAMTMP"/idempotent.bam --sort --log-level FATAL
  $ sorttags "$CRAMTMP"/sorted.bam | sort > "$CRAMTMP"/sorted.sam
  $ sorttags "$CRAMTMP"/idempotent.bam | sort > "$CRAMTMP"/idempotent.sam
  $ diff "$CRAMTMP"/sorted.sam "$CRAMTMP"/idempotent.sam
