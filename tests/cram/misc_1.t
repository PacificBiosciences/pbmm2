  $ IN="$TESTDIR"/data/median.bam
  $ samtools view -h "$IN" > "$CRAMTMP"/sub.sam
  $ samtools view "$IN" | head -n 1 >> "$CRAMTMP"/sub.sam
  $ samtools view -bS "$CRAMTMP"/sub.sam > "$CRAMTMP"/median.bam
  $ IN="$CRAMTMP"/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$REF" "$IN" "$CRAMTMP"/global.alignmentset.xml --log-level FATAL --preset SUBREAD
  $ grep PacBio.AlignmentFile.AlignmentBamFile "$CRAMTMP"/global.alignmentset.xml
  */global.bam* (glob)

  $ "$PBMM2" align -j 1 "$REF" "$IN" local.alignmentset.xml --log-level FATAL --preset SUBREAD
  $ grep PacBio.AlignmentFile.AlignmentBamFile local.alignmentset.xml
  */local.bam* (glob)

  $ mkdir -p "$CRAMTMP"/sub/dir/foo
  $ "$PBMM2" align -j 1 "$REF" "$IN" "$CRAMTMP"/sub/dir/foo/test.alignmentset.xml --log-level FATAL --preset SUBREAD
  $ grep PacBio.AlignmentFile.AlignmentBamFile "$CRAMTMP"/sub/dir/foo/test.alignmentset.xml
  */sub/dir/foo/test.bam* (glob)

  $ mkdir -p "$CRAMTMP"/bla/other/bar
  $ cd "$CRAMTMP"/sub/dir/foo
  $ "$PBMM2" align -j 1 "$REF" "$IN" ../../../bla/other/bar/test.alignmentset.xml --log-level FATAL --preset SUBREAD
  $ grep PacBio.AlignmentFile.AlignmentBamFile "$CRAMTMP"/bla/other/bar/test.alignmentset.xml
  */bla/other/bar/test.bam* (glob)
