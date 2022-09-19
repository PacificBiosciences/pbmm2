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

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/warn_bai.bam --no-bai
  *Option --no-bai has no effect without option --sort!* (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -o -2
  *Gap options have to be strictly positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -O -2
  *Gap options have to be strictly positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -e -2
  *Gap options have to be strictly positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -E -2
  *Gap options have to be strictly positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -k 0
  *Index parameter -k and -w must be positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -k -2
  *Index parameter -k and -w must be positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -w 0
  *Index parameter -k and -w must be positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -w -2
  *Index parameter -k and -w must be positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -o 150
  *Violation of dual gap penalties, E1>E2 and O1+E1<O2+E2* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -z 100 -Z 200
  *Z-drop should not be less than inversion-Z-drop* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -o 20 -O 5 -e 1 -E 2
  *Violation of dual gap penalties, E1>E2 and O1+E1<O2+E2* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam --best-n -1
  *Parameter --best-n, -N must be positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -N -2
  *Parameter --best-n, -N must be positive.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$REF" $BAM $BAM "$CRAMTMP"/fail.bam
  *Incorrect number of arguments. Accepted are at most three!* (glob)
  [1]

  $ FOFN="$CRAMTMP"/fa_not_exist.fofn
  $ echo "FastaNotExist.fasta" > $FOFN
  $ "$PBMM2" align -j 1 "$REF" $FOFN  "$CRAMTMP"/fa_not_exist.bam --preset CCS --rg '@RG\tID:myid\tSM:mysample'
  *Input fofn contains non-existing file: FastaNotExist.fasta* (glob)
  [1]
