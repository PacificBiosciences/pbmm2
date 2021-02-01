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

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ecoli.mmi
  $ REF="$CRAMTMP"/ecoli.mmi

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/read_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$REF" "$IN" "$CRAMTMP"/ref_read.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$REF" "$TESTDIR"/data/median.subreadset.xml "$CRAMTMP"/ref_xml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.subreadset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$TESTDIR"/data/median.subreadset.xml "$REF" "$CRAMTMP"/xml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.subreadset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$REF" "$TESTDIR"/data/median.transcriptset.xml "$CRAMTMP"/ref_transxml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.transcriptset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$TESTDIR"/data/median.transcriptset.xml "$REF" "$CRAMTMP"/transxml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.transcriptset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$REF" "$TESTDIR"/data/median.consensusreadset.xml "$CRAMTMP"/ref_ccsxml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.consensusreadset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$TESTDIR"/data/median.consensusreadset.xml "$REF" "$CRAMTMP"/ccsxml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.consensusreadset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/ref_fasta.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fasta* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/ref_fastq.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fastq* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ "$PBMM2" align -j 1 "$FASTQ" "$REF" "$CRAMTMP"/fastq_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fastq* (glob)
  *REF  input file: *ecoli.mmi* (glob)
