  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/ccs_parameters.bam --log-level DEBUG --preset HiFi 2>&1| grep DEBUG
  *Minimap2 parameters* (glob)
  *Kmer size              : 19 (glob)
  *Minimizer window size  : 10 (glob)
  *Homopolymer compressed : false (glob)
  *Gap open 1             : 5 (glob)
  *Gap open 2             : 56 (glob)
  *Gap extension 1        : 4 (glob)
  *Gap extension 2        : 1 (glob)
  *Match score            : 2 (glob)
  *Mismatch penalty       : 5 (glob)
  *Z-drop                 : 400 (glob)
  *Z-drop inv             : 50 (glob)
  *Bandwidth              : 2000 (glob)
  *Max gap                : 5000 (glob)
  *Long join flank ratio  : 0.5 (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/ccs_parameters.bam --preset foo 2>&1
  *Could not find --preset foo* (glob)
  [1]

Test bam_sort
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/sorted_small.bam --preset SUBREAD --sort -J 1 -m 1M --log-level INFO --log-file "$CRAMTMP"/sorted_small.txt
