  $ IN=$TESTDIR/data/median.bam
  $ samtools view -h $IN > $CRAMTMP/sub.sam
  $ samtools view $IN | head -n 1 >> $CRAMTMP/sub.sam
  $ samtools view -bS $CRAMTMP/sub.sam > $CRAMTMP/median.bam
  $ IN=$CRAMTMP/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ BAM=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ samtools view ${BAM} | awk '{ print "@"$1"\n"$10"\n+\n"$11 }' > $CRAMTMP/median.fastq
  $ FASTQ=$CRAMTMP/median.fastq
  $ samtools view ${BAM} | awk '{ print ">"$1"\n"$10 }' > $CRAMTMP/median.fasta
  $ FASTA=$CRAMTMP/median.fasta

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
    E. Align CCS fastq input and sort output
      $ pbmm2 align ref.fasta movie.Q20.fastq ref.movie.bam --preset CCS --sort --rg '@RG\tID:myid\tSM:mysample'

  $ $__PBTEST_PBMM2_EXE --help 2>&1 | head -n 1
  pbmm2 - minimap2 with native PacBio BAM support* (glob)

  $ $__PBTEST_PBMM2_EXE bla
  ERROR: Unknown tool bla* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align 2>&1 | grep Usage
  Usage: pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml|fa|fq> [out.aligned.bam|xml]* (glob)

  $ $__PBTEST_PBMM2_EXE align --help 2>&1 | grep Usage
  Usage: pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml|fa|fq> [out.aligned.bam|xml]* (glob)

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
  *Could not determine read input type(s). Please do not mix data types, such as BAM+FASTQ. File of files may only contain BAMs or datasets.* (glob)
  [1]

  $ echo $BAM > $CRAMTMP/mixed-bam-fa.fofn
  $ echo $FASTA >> $CRAMTMP/mixed-bam-fa.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-bam-fa.fofn $REF $CRAMTMP/mixed-bam-fq.bam
  *Could not determine read input type(s). Please do not mix data types, such as BAM+FASTQ. File of files may only contain BAMs or datasets.* (glob)
  [1]

  $ echo $FASTQ > $CRAMTMP/mixed-fq-fa.fofn
  $ echo $FASTA >> $CRAMTMP/mixed-fq-fa.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fq-fa.fofn $REF $CRAMTMP/mixed-fq-fa.bam
  *Could not determine read input type(s). Please do not mix data types, such as BAM+FASTQ. File of files may only contain BAMs or datasets.* (glob)
  [1]

  $ echo $FASTQ > $CRAMTMP/mixed-fq-fq.fofn
  $ echo $FASTQ >> $CRAMTMP/mixed-fq-fq.fofn
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/mixed-fq-fq.fofn $REF $CRAMTMP/mixed-fq-fq.bam
  *Could not determine read input type(s). Please do not mix data types, such as BAM+FASTQ. File of files may only contain BAMs or datasets.* (glob)
  [1]
