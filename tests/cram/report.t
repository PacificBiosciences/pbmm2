  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 --report-json mapping_stats.report.json "$IN" "$REF" "$CRAMTMP"/unsorted.bam --preset SUBREAD
  $ grep -c mapped_reads_n mapping_stats.report.json
  1
