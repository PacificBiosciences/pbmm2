  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/default_parameters.bam --log-level DEBUG --preset SUBREAD 2>&1| grep DEBUG
  *Minimap2 parameters* (glob)
  *Kmer size              : 19 (glob)
  *Minimizer window size  : 10 (glob)
  *Homopolymer compressed : true (glob)
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
