  $ IN=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN 2>&1
  *Please provide at least the input arguments: input reference output!* (glob)
  *EXAMPLE: pbmm2 input.subreads.bam reference.fasta output.bam* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN.bam $REF $CRAMTMP/fail.bam 2>&1; rm -rf $CRAMTMP/fail.bam
  *Input data file does not exist* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF.fasta $CRAMTMP/fail.bam 2>&1; rm -rf $CRAMTMP/fail.bam
  *file does not exist* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bumms 2>&1
  *Unknown file extension for output* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $IN $CRAMTMP/fail.bam 2>&1; rm -rf $CRAMTMP/fail.bam
  *Both input files are of type READ. Please check your inputs.* (glob)

  $ $__PBTEST_PBMM2_EXE align $REF $REF $CRAMTMP/fail.bam 2>&1; rm -rf $CRAMTMP/fail.bam
  *Unsupported input data file* (glob)

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

  $ $__PBTEST_PBMM2_EXE --help 2>&1 | head -n 1
  pbmm2 - minimap2 with native PacBio BAM support* (glob)

  $ $__PBTEST_PBMM2_EXE bla
  ERROR: Unknown tool bla* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align 2>&1 | grep Usage
  Usage: pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml> [out.aligned.bam|xml]* (glob)

  $ $__PBTEST_PBMM2_EXE align --help 2>&1 | grep Usage
  Usage: pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml> [out.aligned.bam|xml]* (glob)

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
