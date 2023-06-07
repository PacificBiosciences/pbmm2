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
  *Long join flank ratio  : 0.5 (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/default_overrides.bam --log-level DEBUG --preset SUBREAD -o 5 -O 56 -e 4 -E 1 -k 19 -w 10 -A 2 -B 5 -z 400 -Z 50 -r 1000 -L 0.4 -g 10000 2>&1| grep DEBUG
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

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/ccs_parameters.bam --log-level DEBUG --preset CCS 2>&1| grep DEBUG
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

Test that median filter does not fail
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/median_output.bam --preset SUBREAD --median-filter

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/bestn1.bam --preset SUBREAD --best-n 1
  $ samtools view "$CRAMTMP"/bestn1.bam | wc -l | tr -d ' '
  52

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/rle.bam --preset SUBREAD --collapse-homopolymers
  $ samtools view -H "$CRAMTMP"/rle.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/rle.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/rle.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/rle.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/rle.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/rle.json 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/rle.ref.collapsed.fasta 2> /dev/null | wc -l | tr -d ' '
  1
