  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align "$IN" "$REF" "$CRAMTMP"/default_overrides.bam --log-level DEBUG -o 5 -O 56 -e 4 -E 1 -k 19 -w 10 -A 2 -B 5 -z 400 -Z 50 -r 1000 -L 0.4 -g 10000 2>&1| grep DEBUG
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
  *Bandwidth              : 1000 (glob)
  *Max gap                : 10000 (glob)
  *Long join flank ratio  : 0.4 (glob)
