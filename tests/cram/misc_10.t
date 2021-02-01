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

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset CCS --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: CCS / HiFi (glob)
  * Kmer size              : 19 (glob)
  * Minimizer window size  : 10 (glob)
  * Homopolymer compressed : false (glob)

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset hifi --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: CCS / HiFi (glob)
  * Kmer size              : 19 (glob)
  * Minimizer window size  : 10 (glob)
  * Homopolymer compressed : false (glob)

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset HiFi --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: CCS / HiFi (glob)
  * Kmer size              : 19 (glob)
  * Minimizer window size  : 10 (glob)
  * Homopolymer compressed : false (glob)

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset isoseq --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: ISOSEQ (glob)
  * Kmer size              : 15 (glob)
  * Minimizer window size  : 5 (glob)
  * Homopolymer compressed : false (glob)

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset sUBreAd --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: SUBREAD (glob)
  * Kmer size              : 19 (glob)
  * Minimizer window size  : 10 (glob)
  * Homopolymer compressed : true (glob)
