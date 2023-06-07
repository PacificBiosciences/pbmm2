  $ IN="$TESTDIR"/data/m54019_171011_032401_tiny.subreadset.xml
  $ REF="$TESTDIR"/data/lambdaNEB_BsaAI_allFrags_incLeftRightEnds_unrolled_250k.fasta
  $ source "$TESTDIR"/sorttags.sh

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/sorted.bam --preset SUBREAD --sort --log-level FATAL --zmw
  $ "$PBMM2" align -j 1 "$CRAMTMP"/sorted.bam "$REF" "$CRAMTMP"/idempotent.bam --preset SUBREAD --sort --zmw --log-level FATAL
  $ sorttags "$CRAMTMP"/sorted.bam | sort > "$CRAMTMP"/sorted.sam
  $ sorttags "$CRAMTMP"/idempotent.bam | sort > "$CRAMTMP"/idempotent.sam
  $ diff "$CRAMTMP"/sorted.sam "$CRAMTMP"/idempotent.sam
