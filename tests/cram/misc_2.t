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

  $ "$PBMM2" align -j 1 "$IN" 2>&1
  *Please provide at least the input arguments: reference input output!* (glob)
  *EXAMPLE: pbmm2 reference.fasta input.subreads.bam output.bam* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN".bam "$REF" "$CRAMTMP"/fail.bam 2>&1; rm -rf "$CRAMTMP"/fail.bam
  *Input data file does not exist* (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF".fasta "$CRAMTMP"/fail.bam 2>&1; rm -rf "$CRAMTMP"/fail.bam
  *file does not exist* (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bumms 2>&1
  *Unknown file extension for output* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$IN" "$CRAMTMP"/fail.bam 2>&1; rm -rf "$CRAMTMP"/fail.bam
  *Both input files are of type BAM. Please check your inputs.* (glob)

  $ "$PBMM2" align -j 1 "$REF" "$REF" "$CRAMTMP"/fail.bam 2>&1; rm -rf "$CRAMTMP"/fail.bam
  *Input is FASTA.* (glob)

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam -L 1.1
  *Option -L,--lj-min-ratio has to be between a ratio betweem 0 and 1.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam --zmw
  *Option --zmw can only be used with a subreadset.xml containing subread + scraps BAM files.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam --hqregion
  *Option --hqregion can only be used with a subreadset.xml containing subread + scraps BAM files.* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/fail.bam --zmw --hqregion 2>&1; rm -rf "$CRAMTMP"/fail.bam
  *Options --zmw, --hqregion and --median-filter are mutually exclusive.* (glob)

  $ "$PBMM2" align "$IN" "$REF" "$CRAMTMP"/fail.bam --sort -J 1 -m 1000P 2>&1; rm -rf "$CRAMTMP"/fail.bam
  *Unknown size multiplier P* (glob)

  $ "$PBMM2" align "$IN" "$REF" "$CRAMTMP"/fail.bam --sort -J 1 -m 10000000000 2>&1; rm -rf "$CRAMTMP"/fail.bam

  $ "$PBMM2" index "$REF" "$CRAMTMP"/index_logging.mmi --log-file "$CRAMTMP"/index_logging.txt 2>&1

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ref.mmi
  $ "$PBMM2" align "$CRAMTMP"/ref.mmi "$IN" "$CRAMTMP"/mmi.fail.bam --collapse-homopolymers
  *Cannot combine --collapse-homopolymers with MMI input.* (glob)
  [1]

  $ "$PBMM2" index "$REF" 2>&1
  *Please provide both arguments: input output!* (glob)
  *EXAMPLE: pbmm2 index reference.fasta output.mmi* (glob)
  [1]

  $ "$PBMM2" index "$REF" "$CRAMTMP"/fail.mmi "$CRAMTMP"/fail.mmi 2>&1
  *Please provide both arguments: input output!* (glob)
  *EXAMPLE: pbmm2 index reference.fasta output.mmi* (glob)
  [1]

  $ "$PBMM2" index "$REF".fasta 2>&1
  *Please provide both arguments: input output!* (glob)
  *EXAMPLE: pbmm2 index reference.fasta output.mmi* (glob)
  [1]

  $ "$PBMM2" index "$IN" "$CRAMTMP"/fail.mmi 2>&1
  *Unsupported input data file* (glob)
  [1]

  $ "$PBMM2" index "$REF" "$CRAMTMP"/fail.mmx 2>&1
  *Output file must end with .mmi:* (glob)
  [1]

  $ "$PBMM2"
  pbmm2 - minimap2 with native PacBio BAM support
  * (glob)
  Usage:
    pbmm2 <tool>
  * (glob)
    -h,--help    Show this help and exit.
    --version    Show application version and exit.
  * (glob)
  Tools:
    index      Index reference and store as .mmi file
    align      Align PacBio reads to reference sequences
  * (glob)
  Examples:
    pbmm2 index ref.referenceset.xml ref.mmi
    pbmm2 align ref.referenceset.xml movie.subreadset.xml ref.movie.alignmentset.xml
  * (glob)
  Typical workflows:
    A. Generate index file for reference and reuse it to align reads
       $ pbmm2 index ref.fasta ref.mmi
       $ pbmm2 align ref.mmi movie.subreads.bam ref.movie.bam
  * (glob)
    B. Align reads and sort on-the-fly, with 4 alignment and 2 sort threads
       $ pbmm2 align ref.fasta movie.subreads.bam ref.movie.bam --sort -j 4 -J 2
  * (glob)
    C. Align reads, sort on-the-fly, and create PBI
       $ pbmm2 align ref.fasta movie.subreadset.xml ref.movie.alignmentset.xml --sort
  * (glob)
    D. Omit output file and stream BAM output to stdout
       $ pbmm2 align hg38.mmi movie1.subreadset.xml | samtools sort > hg38.movie1.sorted.bam
  * (glob)
    E. Align CCS fastq input and sort on-the-fly
       $ pbmm2 align ref.fasta movie.Q20.fastq ref.movie.bam --preset CCS --sort --rg '@RG\tID:myid\tSM:mysample'

  $ "$PBMM2" --help 2>&1 | head -n 1
  pbmm2 - minimap2 with native PacBio BAM support* (glob)

  $ "$PBMM2" bla
  pbmm2 ERROR: [pbcopper] command line ERROR: unknown tool 'bla' requested (glob)
  [1]

  $ "$PBMM2" align 2>&1 | grep "\[options\]"
    pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml|fa|fq|gz|fofn> [out.aligned.bam|xml]* (glob)

  $ "$PBMM2" align --help 2>&1 | grep "\[options\]"
    pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml|fa|fq|gz|fofn> [out.aligned.bam|xml]* (glob)

  $ "$PBMM2" index 2>&1 | grep "\[options\]"
    pbmm2 index [options] <ref.fa|xml> <out.mmi>* (glob)

  $ "$PBMM2" --version
  pbmm2 *.*.* (*) (glob)
  $ "$PBMM2" align --version
  pbmm2 *.*.* (*) (glob)
  $ "$PBMM2" index --version
  pbmm2 *.*.* (*) (glob)
