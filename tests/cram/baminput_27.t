  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/ccs_parameters.bam --log-level DEBUG --preset HiFi 2>&1| grep DEBUG
  *Minimap2 parameters* (glob)
  *Kmer size              : 19 (glob)
  *Minimizer window size  : 19 (glob)
  *Homopolymer compressed : false (glob)
  *Gap open 1             : 6 (glob)
  *Gap open 2             : 26 (glob)
  *Gap extension 1        : 2 (glob)
  *Gap extension 2        : 1 (glob)
  *Match score            : 1 (glob)
  *Mismatch penalty       : 4 (glob)
  *Z-drop                 : 400 (glob)
  *Z-drop inv             : 50 (glob)
  *Bandwidth              : 2000 (glob)
  *Max gap                : 10000 (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/ccs_parameters.bam --preset foo 2>&1
  *Could not find --preset foo* (glob)
  [1]

Test bam_sort
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/sorted_small.bam --preset SUBREAD --sort -J 1 -m 1M --log-level INFO --log-file "$CRAMTMP"/sorted_small.txt
