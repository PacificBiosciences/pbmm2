  $ IN="$TESTDIR"/data/median.bam
  $ "$SAMTOOLS" view -h "$IN" > "$CRAMTMP"/sub.sam
  $ "$SAMTOOLS" view "$IN" | head -n 1 >> "$CRAMTMP"/sub.sam
  $ "$SAMTOOLS" view -bS "$CRAMTMP"/sub.sam > "$CRAMTMP"/median.bam
  $ IN="$CRAMTMP"/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ BAM="$TESTDIR"/data/median.bam
  $ "$SAMTOOLS" view "$BAM" | awk '{ print "@"$1"\n"$10"\n+\n"$11 }' > "$CRAMTMP"/median.fastq
  $ FASTQ="$CRAMTMP"/median.fastq
  $ cp "$CRAMTMP"/median.fastq "$CRAMTMP"/median_compressed.fastq
  $ gzip "$CRAMTMP"/median_compressed.fastq
  $ FASTQGZ="$CRAMTMP"/median_compressed.fastq.gz
  $ "$SAMTOOLS" view "$BAM" | awk '{ print ">"$1"\n"$10 }' > "$CRAMTMP"/median.fasta
  $ FASTA="$CRAMTMP"/median.fasta
  $ cp "$CRAMTMP"/median.fasta "$CRAMTMP"/median_compressed.fasta
  $ gzip "$CRAMTMP"/median_compressed.fasta
  $ FASTAGZ="$CRAMTMP"/median_compressed.fasta.gz

  $ "$PBMM2" align -j 1 "$REF" "$IN" "$CRAMTMP"/global.alignmentset.xml --log-level FATAL
  $ grep PacBio.AlignmentFile.AlignmentBamFile "$CRAMTMP"/global.alignmentset.xml
  */global.bam* (glob)

  $ "$PBMM2" align -j 1 "$REF" "$IN" local.alignmentset.xml --log-level FATAL
  $ grep PacBio.AlignmentFile.AlignmentBamFile local.alignmentset.xml
  */local.bam* (glob)

  $ mkdir -p "$CRAMTMP"/sub/dir/foo
  $ "$PBMM2" align -j 1 "$REF" "$IN" "$CRAMTMP"/sub/dir/foo/test.alignmentset.xml --log-level FATAL
  $ grep PacBio.AlignmentFile.AlignmentBamFile "$CRAMTMP"/sub/dir/foo/test.alignmentset.xml
  */sub/dir/foo/test.bam* (glob)

  $ mkdir -p "$CRAMTMP"/bla/other/bar
  $ cd "$CRAMTMP"/sub/dir/foo
  $ "$PBMM2" align -j 1 "$REF" "$IN" ../../../bla/other/bar/test.alignmentset.xml --log-level FATAL
  $ grep PacBio.AlignmentFile.AlignmentBamFile "$CRAMTMP"/bla/other/bar/test.alignmentset.xml
  */bla/other/bar/test.bam* (glob)

  $ cd ../../../
