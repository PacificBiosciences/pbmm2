  $ BAM="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

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

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/fastq_unsorted.bam 2>&1
  *Input is FASTQ.* (glob)
  $ samtools view -H "$CRAMTMP"/fastq_unsorted.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsorted.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsorted.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsorted.json 2> /dev/null | wc -l | tr -d ' '
  0
