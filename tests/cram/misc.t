  $ IN=$TESTDIR/data/median.bam
  $ samtools view -h $IN > $CRAMTMP/sub.sam
  $ samtools view $IN | head -n 1 >> $CRAMTMP/sub.sam
  $ samtools view -bS $CRAMTMP/sub.sam > $CRAMTMP/median.bam
  $ IN=$CRAMTMP/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ BAM=$TESTDIR/data/median.bam
  $ samtools view ${BAM} | awk '{ print "@"$1"\n"$10"\n+\n"$11 }' > $CRAMTMP/median.fastq
  $ FASTQ=$CRAMTMP/median.fastq
  $ cp $CRAMTMP/median.fastq $CRAMTMP/median_compressed.fastq
  $ gzip $CRAMTMP/median_compressed.fastq
  $ FASTQGZ=$CRAMTMP/median_compressed.fastq.gz
  $ samtools view ${BAM} | awk '{ print ">"$1"\n"$10 }' > $CRAMTMP/median.fasta
  $ FASTA=$CRAMTMP/median.fasta
  $ cp $CRAMTMP/median.fasta $CRAMTMP/median_compressed.fasta
  $ gzip $CRAMTMP/median_compressed.fasta
  $ FASTAGZ=$CRAMTMP/median_compressed.fasta.gz

  $ $__PBTEST_PBMM2_EXE align $REF $IN $CRAMTMP/global.alignmentset.xml --log-level FATAL
  $ grep PacBio.AlignmentFile.AlignmentBamFile $CRAMTMP/global.alignmentset.xml
  */global.bam* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $IN local.alignmentset.xml --log-level FATAL
  $ grep PacBio.AlignmentFile.AlignmentBamFile local.alignmentset.xml
  */local.bam* (glob)

  $ mkdir -p $CRAMTMP/sub/dir/foo
  $ $__PBTEST_PBMM2_EXE align $REF $IN $CRAMTMP/sub/dir/foo/test.alignmentset.xml --log-level FATAL
  $ grep PacBio.AlignmentFile.AlignmentBamFile $CRAMTMP/sub/dir/foo/test.alignmentset.xml
  */sub/dir/foo/test.bam* (glob)

  $ mkdir -p $CRAMTMP/bla/other/bar
  $ cd $CRAMTMP/sub/dir/foo
  $ $__PBTEST_PBMM2_EXE align $REF $IN ../../../bla/other/bar/test.alignmentset.xml --log-level FATAL
  $ grep PacBio.AlignmentFile.AlignmentBamFile $CRAMTMP/bla/other/bar/test.alignmentset.xml
  */bla/other/bar/test.bam* (glob)

  $ cd ../../../

  $ $__PBTEST_PBMM2_EXE align $IN 2>&1
  *Please provide at least the input arguments: reference input output!* (glob)
  *EXAMPLE: pbmm2 reference.fasta input.subreads.bam output.bam* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN.bam $REF $CRAMTMP/fail.bam 2>&1; rm -rf $CRAMTMP/fail.bam
  *Input data file does not exist* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF.fasta $CRAMTMP/fail.bam 2>&1; rm -rf $CRAMTMP/fail.bam
  *file does not exist* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bumms 2>&1
  *Unknown file extension for output* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $IN $CRAMTMP/fail.bam 2>&1; rm -rf $CRAMTMP/fail.bam
  *Both input files are of type BAM. Please check your inputs.* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $REF $CRAMTMP/fail.bam 2>&1; rm -rf $CRAMTMP/fail.bam
  *Input is FASTA.* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -L 1.1
  *Option -L,--lj-min-ratio has to be between a ratio betweem 0 and 1.* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam --zmw
  *Option --zmw can only be used with a subreadset.xml containing subread + scraps BAM files.* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam --hqregion
  *Option --hqregion can only be used with a subreadset.xml containing subread + scraps BAM files.* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam --zmw --hqregion 2>&1; rm -rf $CRAMTMP/fail.bam
  *Options --zmw, --hqregion and --median-filter are mutually exclusive.* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam --sort -J 1 -m 1000G 2>&1; rm -rf $CRAMTMP/fail.bam
  *Trying to allocate more memory for sorting* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam --sort -J 1 -m 1000P 2>&1; rm -rf $CRAMTMP/fail.bam
  *Unknown size multiplier P* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam --sort -J 1 -m 10000000000 2>&1; rm -rf $CRAMTMP/fail.bam

  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/index_logging.mmi --log-file $CRAMTMP/index_logging.txt 2>&1

  $ $__PBTEST_PBMM2_EXE index $REF 2>&1
  *Please provide both arguments: input output!* (glob)
  *EXAMPLE: pbmm2 index reference.fasta output.mmi* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/fail.mmi $CRAMTMP/fail.mmi 2>&1
  *Please provide both arguments: input output!* (glob)
  *EXAMPLE: pbmm2 index reference.fasta output.mmi* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE index $REF.fasta 2>&1
  *Please provide both arguments: input output!* (glob)
  *EXAMPLE: pbmm2 index reference.fasta output.mmi* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE index $IN $CRAMTMP/fail.mmi 2>&1
  *Unsupported input data file* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/fail.mmx 2>&1
  *Output file must end with .mmi:* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE
  pbmm2 - minimap2 with native PacBio BAM support
  * (glob)
  Usage:
    pbmm2 <tool>
  * (glob)
  Options:
    -h, --help   Output this help.
    --version    Output version info.
  * (glob)
  Tools:
      index      Index reference and store as .mmi file
      align      Align PacBio reads to reference sequences
  * (glob)
  Examples:
    pbmm2 align ref.referenceset.xml movie.subreadset.xml ref.movie.alignmentset.xml
    pbmm2 index ref.referenceset.xml ref.mmi
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

  $ $__PBTEST_PBMM2_EXE --help 2>&1 | head -n 1
  pbmm2 - minimap2 with native PacBio BAM support* (glob)

  $ $__PBTEST_PBMM2_EXE bla
  ERROR: Unknown tool bla* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align 2>&1 | grep Usage
  Usage: pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml|fa|fq|gz|fofn> [out.aligned.bam|xml]* (glob)

  $ $__PBTEST_PBMM2_EXE align --help 2>&1 | grep Usage
  Usage: pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml|fa|fq|gz|fofn> [out.aligned.bam|xml]* (glob)

  $ $__PBTEST_PBMM2_EXE index 2>&1 | grep Usage
  Usage: pbmm2 index [options] <ref.fa|xml> <out.mmi>* (glob)

  $ $__PBTEST_PBMM2_EXE --version
  pbmm2 *.*.* (*) (glob)
  $ $__PBTEST_PBMM2_EXE align --version
  pbmm2 *.*.* (*) (glob)
  $ $__PBTEST_PBMM2_EXE index --version
  pbmm2 *.*.* (*) (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/read_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $IN $CRAMTMP/ref_read.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $TESTDIR/data/median.subreadset.xml $CRAMTMP/ref_xml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.subreadset.xml* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.subreadset.xml $REF $CRAMTMP/xml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.subreadset.xml* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $TESTDIR/data/median.transcriptset.xml $CRAMTMP/ref_transxml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.transcriptset.xml* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.transcriptset.xml $REF $CRAMTMP/transxml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.transcriptset.xml* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $TESTDIR/data/median.consensusreadset.xml $CRAMTMP/ref_ccsxml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.consensusreadset.xml* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.consensusreadset.xml $REF $CRAMTMP/ccsxml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.consensusreadset.xml* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/ref_fasta.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fasta* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTQ $CRAMTMP/ref_fastq.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fastq* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE align $FASTQ $REF $CRAMTMP/fastq_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fastq* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)

  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/ecoli.mmi
  $ REF=$CRAMTMP/ecoli.mmi

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/read_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $IN $CRAMTMP/ref_read.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $TESTDIR/data/median.subreadset.xml $CRAMTMP/ref_xml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.subreadset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.subreadset.xml $REF $CRAMTMP/xml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.subreadset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $TESTDIR/data/median.transcriptset.xml $CRAMTMP/ref_transxml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.transcriptset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.transcriptset.xml $REF $CRAMTMP/transxml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.transcriptset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $TESTDIR/data/median.consensusreadset.xml $CRAMTMP/ref_ccsxml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.consensusreadset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.consensusreadset.xml $REF $CRAMTMP/ccsxml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.consensusreadset.xml* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/ref_fasta.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fasta* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTQ $CRAMTMP/ref_fastq.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fastq* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ $__PBTEST_PBMM2_EXE align $FASTQ $REF $CRAMTMP/fastq_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fastq* (glob)
  *REF  input file: *ecoli.mmi* (glob)

  $ REF=$TESTDIR/data/ecoli.referenceset.xml

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/read_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $IN $CRAMTMP/ref_read.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $TESTDIR/data/median.subreadset.xml $CRAMTMP/ref_xml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.subreadset.xml* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.subreadset.xml $REF $CRAMTMP/xml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.subreadset.xml* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $TESTDIR/data/median.transcriptset.xml $CRAMTMP/ref_transxml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.transcriptset.xml* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.transcriptset.xml $REF $CRAMTMP/transxml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.transcriptset.xml* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $TESTDIR/data/median.consensusreadset.xml $CRAMTMP/ref_ccsxml.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.consensusreadset.xml* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.consensusreadset.xml $REF $CRAMTMP/ccsxml_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.consensusreadset.xml* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/ref_fasta.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fasta* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTQ $CRAMTMP/ref_fastq.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fastq* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ $__PBTEST_PBMM2_EXE align $FASTQ $REF $CRAMTMP/fastq_ref.bam --log-level INFO 2>&1 | grep "input file"
  *READ input file: *median.fastq* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)

  $ echo $BAM > $CRAMTMP/mixed-bam-fq.fofn
  $ echo $FASTQ >> $CRAMTMP/mixed-bam-fq.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-bam-fq.fofn $REF $CRAMTMP/mixed-bam-fq.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]

  $ echo $BAM > $CRAMTMP/mixed-bam-fa.fofn
  $ echo $FASTA >> $CRAMTMP/mixed-bam-fa.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-bam-fa.fofn $REF $CRAMTMP/mixed-bam-fq.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]

  $ echo $FASTQ > $CRAMTMP/mixed-fq-fa.fofn
  $ echo $FASTA >> $CRAMTMP/mixed-fq-fa.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fq-fa.fofn $REF $CRAMTMP/mixed-fq-fa.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]

  $ echo $FASTQGZ > $CRAMTMP/mixed-fqgz-fa.fofn
  $ echo $FASTA >> $CRAMTMP/mixed-fqgz-fa.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fqgz-fa.fofn $REF $CRAMTMP/mixed-fq-fa.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]

  $ echo $FASTQ > $CRAMTMP/mixed-fq-fagz.fofn
  $ echo $FASTAGZ >> $CRAMTMP/mixed-fq-fagz.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fq-fagz.fofn $REF $CRAMTMP/mixed-fq-fa.bam
  *Input fofn contains different file types. This is not supported.* (glob)
  [1]

  $ echo $FASTQ > $CRAMTMP/mixed-fq-fq.fofn
  $ echo $FASTQ >> $CRAMTMP/mixed-fq-fq.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fq-fq.fofn $REF $CRAMTMP/mixed-fq-fq.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTQ FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fq-fq.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Mapped Concordance: 91* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo $FASTA > $CRAMTMP/mixed-fa-fa.fofn
  $ echo $FASTA >> $CRAMTMP/mixed-fa-fa.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fa-fa.fofn $REF $CRAMTMP/mixed-fa-fa.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTA FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fa-fa.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Mapped Concordance: 91* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo $FASTQGZ > $CRAMTMP/mixed-fqgz-fqgz.fofn
  $ echo $FASTQGZ >> $CRAMTMP/mixed-fqgz-fqgz.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fqgz-fqgz.fofn $REF $CRAMTMP/mixed-fqgz-fqgz.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTQ FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fqgz-fqgz.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Mapped Concordance: 91* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo $FASTAGZ > $CRAMTMP/mixed-fagz-fagz.fofn
  $ echo $FASTAGZ >> $CRAMTMP/mixed-fagz-fagz.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fagz-fagz.fofn $REF $CRAMTMP/mixed-fagz-fagz.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTA FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fagz-fagz.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Mapped Concordance: 91* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo $FASTQGZ > $CRAMTMP/mixed-fqgz-fq.fofn
  $ echo $FASTQ >> $CRAMTMP/mixed-fqgz-fq.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fqgz-fq.fofn $REF $CRAMTMP/mixed-fqgz-fq.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTQ FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fqgz-fq.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Mapped Concordance: 91* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo $FASTA > $CRAMTMP/mixed-fa-fagz.fofn
  $ echo $FASTAGZ >> $CRAMTMP/mixed-fa-fagz.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fa-fagz.fofn $REF $CRAMTMP/mixed-fa-fagz.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *Input is FASTA FOFN. Output BAM file cannot be used for polishing with GenomicConsensus!* (glob)
  *READ input file: *mixed-fa-fagz.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Mapped Concordance: 91* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ echo $BAM > $CRAMTMP/mixed-bam-bam.fofn
  $ echo $BAM >> $CRAMTMP/mixed-bam-bam.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-bam-bam.fofn $REF $CRAMTMP/mixed-bam-bam.bam -j 4 --log-level INFO 2>&1
  *Using 4 threads for alignments.* (glob)
  *READ input file: *mixed-bam-bam.fofn* (glob)
  *REF  input file: *ecoli.referenceset.xml* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Mapped Reads: 104 (glob)
  *Alignments: 192 (glob)
  *Mapped Bases: 484874 (glob)
  *Mean Mapped Concordance: 91* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF --sort --no-bai > $CRAMTMP/warn_bai_pipe.bam
  *Option --no-bai has no effect when using an output pipe!* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/warn_bai.bam --no-bai
  *Option --no-bai has no effect without option --sort!* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -o -2
  *Gap options have to be strictly positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -O -2
  *Gap options have to be strictly positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -e -2
  *Gap options have to be strictly positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -E -2
  *Gap options have to be strictly positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -k 0
  *Index parameter -k and -w must be positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -k -2
  *Index parameter -k and -w must be positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -w 0
  *Index parameter -k and -w must be positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -w -2
  *Index parameter -k and -w must be positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -o 150
  *Violation of dual gap penalties, E1>E2 and O1+E1<O2+E2 (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -z 100 -Z 200
  *Z-drop should not be less than inversion-Z-drop (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -o 20 -O 5 -e 1 -E 2
  *Violation of dual gap penalties, E1>E2 and O1+E1<O2+E2 (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam --best-n -1
  *Parameter --best-n, -N must be positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -N -2
  *Parameter --best-n, -N must be positive. (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $REF $BAM $BAM $CRAMTMP/fail.bam
  *Incorrect number of arguments. Accepted are at most three!* (glob)
  [1]
