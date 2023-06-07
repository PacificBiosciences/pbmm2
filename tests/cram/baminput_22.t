  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -j 500 --preset SUBREAD 2>&1| grep WARN
  *Requested more threads for alignment (500) than system-wide available* (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -J 500 --preset SUBREAD --sort 2>&1| grep AlignSettings | grep Requested
  *Requested more threads for sorting (500) and alignment (1) than system-wide available* (glob)
  *Requested more threads for sorting (500) than system-wide available* (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/pass2.bam -j 1 -J 500 --preset SUBREAD -m 500G
  *Requested 500 threads for sorting, without specifying --sort. Please check your input. (glob)
  *Requested 500G memory for sorting, without specifying --sort. Please check your input. (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail3.bam -j 1 -J 500 --preset SUBREAD --sort --sort 2>&1| grep Requested
  *Requested more threads for sorting* (glob)
  *Requested more threads for sorting* (glob)

  $ "$PBMM2" align "$IN" "$REF" "$CRAMTMP"/sort_percentage_4.bam -j 4 --preset SUBREAD --sort --log-level INFO 2>&1 | grep "threads for alignments"
  *Using 3 threads for alignments, 1 threads for sorting, and 768M bytes RAM for sorting.* (glob)
