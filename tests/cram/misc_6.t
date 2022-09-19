  $ IN="$TESTDIR"/data/median.bam
  $ samtools view -h "$IN" > "$CRAMTMP"/sub.sam
  $ samtools view "$IN" | head -n 1 >> "$CRAMTMP"/sub.sam
  $ samtools view -bS "$CRAMTMP"/sub.sam > "$CRAMTMP"/median.bam
  $ IN="$CRAMTMP"/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ BAM="$TESTDIR"/data/median.bam
  $ samtools view "$BAM" | awk '{ print "@"$1"\n"$10"\n+\n"$11 }' > "$CRAMTMP"/median.fastq
  $ FASTQ="$CRAMTMP"/median.fastq
  $ cp "$CRAMTMP"/median.fastq "$CRAMTMP"/median_compressed.fastq
  $ gzip "$CRAMTMP"/median_compressed.fastq
  $ FASTQGZ="$CRAMTMP"/median_compressed.fastq.gz
  $ samtools view "$BAM" | awk '{ print ">"$1"\n"$10 }' > "$CRAMTMP"/median.fasta
  $ FASTA="$CRAMTMP"/median.fasta
  $ cp "$CRAMTMP"/median.fasta "$CRAMTMP"/median_compressed.fasta
  $ gzip "$CRAMTMP"/median_compressed.fasta
  $ FASTAGZ="$CRAMTMP"/median_compressed.fasta.gz

  $ REF="$TESTDIR"/data/ecoli.referenceset.xml

  $ echo $BAM > "$CRAMTMP"/mixed-bam-fq.fofn
  $ echo "$FASTQ" >> "$CRAMTMP"/mixed-bam-fq.fofn
  $ "$PBMM2" align -j 1 "$CRAMTMP"/mixed-bam-fq.fofn "$REF" "$CRAMTMP"/mixed-bam-fq.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]

  $ echo $BAM > "$CRAMTMP"/mixed-bam-fa.fofn
  $ echo "$FASTA" >> "$CRAMTMP"/mixed-bam-fa.fofn
  $ "$PBMM2" align -j 1 "$CRAMTMP"/mixed-bam-fa.fofn "$REF" "$CRAMTMP"/mixed-bam-fq.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]

  $ echo "$FASTQ" > "$CRAMTMP"/mixed-fq-fa.fofn
  $ echo "$FASTA" >> "$CRAMTMP"/mixed-fq-fa.fofn
  $ "$PBMM2" align -j 1 "$CRAMTMP"/mixed-fq-fa.fofn "$REF" "$CRAMTMP"/mixed-fq-fa.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]

  $ echo "$FASTQGZ" > "$CRAMTMP"/mixed-fqgz-fa.fofn
  $ echo "$FASTA" >> "$CRAMTMP"/mixed-fqgz-fa.fofn
  $ "$PBMM2" align -j 1 "$CRAMTMP"/mixed-fqgz-fa.fofn "$REF" "$CRAMTMP"/mixed-fq-fa.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]

  $ echo "$FASTQ" > "$CRAMTMP"/mixed-fq-fagz.fofn
  $ echo "$FASTAGZ" >> "$CRAMTMP"/mixed-fq-fagz.fofn
  $ "$PBMM2" align -j 1 "$CRAMTMP"/mixed-fq-fagz.fofn "$REF" "$CRAMTMP"/mixed-fq-fa.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]
