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

  $ echo "$FASTA" > "$CRAMTMP"/mixed-fa-fa.fofn
  $ echo "$FASTA" >> "$CRAMTMP"/mixed-fa-fa.fofn
  $ "$PBMM2" align "$CRAMTMP"/mixed-fa-fa.fofn "$REF" "$CRAMTMP"/mixed-fa-fa.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTA FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fa-fa.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Gap-Compressed Sequence Identity* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo "$FASTQGZ" > "$CRAMTMP"/mixed-fqgz-fqgz.fofn
  $ echo "$FASTQGZ" >> "$CRAMTMP"/mixed-fqgz-fqgz.fofn
  $ "$PBMM2" align "$CRAMTMP"/mixed-fqgz-fqgz.fofn "$REF" "$CRAMTMP"/mixed-fqgz-fqgz.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTQ FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fqgz-fqgz.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Gap-Compressed Sequence Identity* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo "$FASTAGZ" > "$CRAMTMP"/mixed-fagz-fagz.fofn
  $ echo "$FASTAGZ" >> "$CRAMTMP"/mixed-fagz-fagz.fofn
  $ "$PBMM2" align "$CRAMTMP"/mixed-fagz-fagz.fofn "$REF" "$CRAMTMP"/mixed-fagz-fagz.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTA FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fagz-fagz.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Gap-Compressed Sequence Identity* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo "$FASTQGZ" > "$CRAMTMP"/mixed-fqgz-fq.fofn
  $ echo "$FASTQ" >> "$CRAMTMP"/mixed-fqgz-fq.fofn
  $ "$PBMM2" align "$CRAMTMP"/mixed-fqgz-fq.fofn "$REF" "$CRAMTMP"/mixed-fqgz-fq.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTQ FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fqgz-fq.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Gap-Compressed Sequence Identity* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo "$FASTA" > "$CRAMTMP"/mixed-fa-fagz.fofn
  $ echo "$FASTAGZ" >> "$CRAMTMP"/mixed-fa-fagz.fofn
  $ "$PBMM2" align "$CRAMTMP"/mixed-fa-fagz.fofn "$REF" "$CRAMTMP"/mixed-fa-fagz.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTA FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fa-fagz.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Gap-Compressed Sequence Identity* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo $BAM > "$CRAMTMP"/mixed-bam-bam.fofn
  $ echo $BAM >> "$CRAMTMP"/mixed-bam-bam.fofn
  $ "$PBMM2" align "$CRAMTMP"/mixed-bam-bam.fofn "$REF" "$CRAMTMP"/mixed-bam-bam.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *READ input file: *mixed-bam-bam.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Gap-Compressed Sequence Identity* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)
