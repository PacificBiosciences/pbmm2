  $ IN="$TESTDIR"/data/median.bam
  $ samtools view -h "$IN" > "$CRAMTMP"/sub.sam
  $ samtools view "$IN" | head -n 1 >> "$CRAMTMP"/sub.sam
  $ samtools view -bS "$CRAMTMP"/sub.sam > "$CRAMTMP"/median.bam
  $ IN="$CRAMTMP"/median.bam
  $ REF="$TESTDIR"/data/ecoli.referenceset.xml

  $ BAM="$TESTDIR"/data/median.bam

  $ samtools view "$BAM" | awk '{ print "@"$1"\n"$10"\n+\n"$11 }' > "$CRAMTMP"/median.fastq
  $ FASTQ="$CRAMTMP"/median.fastq

  $ echo "$FASTQ" > "$CRAMTMP"/mixed-fq-fq.fofn
  $ echo "$FASTQ" >> "$CRAMTMP"/mixed-fq-fq.fofn
  $ "$PBMM2" align "$CRAMTMP"/mixed-fq-fq.fofn "$REF" "$CRAMTMP"/mixed-fq-fq.bam -j 4 --log-level INFO --preset SUBREAD 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTQ FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fq-fq.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 486732 (glob)
  *Mean Gap-Compressed Sequence Identity* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)
