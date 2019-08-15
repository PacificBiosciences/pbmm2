  $ FASTA=$TESTDIR/data/empty.fasta
  $ FASTQ=$TESTDIR/data/empty.fastq
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta
  $ REFXML=$TESTDIR/data/ecoli.referenceset.xml

FASTA TESTS
  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/empty_fasta_aligned_unsorted.bam -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fasta_aligned_unsorted.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:unknown\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/empty_fasta_aligned_sorted.bam --sort -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fasta_aligned_sorted.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:coordinate\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTA $CRAMTMP/empty_fasta_aligned_sorted_ccs.bam --sort --preset CCS -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fasta_aligned_sorted_ccs.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:coordinate\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REFXML $FASTA $CRAMTMP/empty_fasta_xml_aligned_unsorted.bam -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fasta_xml_aligned_unsorted.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:unknown\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REFXML $FASTA $CRAMTMP/empty_fasta_xml_aligned_sorted.bam --sort -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fasta_xml_aligned_sorted.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:coordinate\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REFXML $FASTA $CRAMTMP/empty_fasta_xml_aligned_sorted_ccs.bam --sort --preset CCS -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fasta_xml_aligned_sorted_ccs.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:coordinate\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

FASTQ tests
  $ $__PBTEST_PBMM2_EXE align $REF $FASTQ $CRAMTMP/empty_fastq_aligned_unsorted.bam -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fastq_aligned_unsorted.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:unknown\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTQ $CRAMTMP/empty_fastq_aligned_sorted.bam --sort -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fastq_aligned_sorted.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:coordinate\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REF $FASTQ $CRAMTMP/empty_fastq_aligned_sorted_ccs.bam --sort --preset CCS -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fastq_aligned_sorted_ccs.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:coordinate\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REFXML $FASTQ $CRAMTMP/empty_fastq_xml_aligned_unsorted.bam -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fastq_xml_aligned_unsorted.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:unknown\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REFXML $FASTQ $CRAMTMP/empty_fastq_xml_aligned_sorted.bam --sort -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fastq_xml_aligned_sorted.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:coordinate\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ $__PBTEST_PBMM2_EXE align $REFXML $FASTQ $CRAMTMP/empty_fastq_xml_aligned_sorted_ccs.bam --sort --preset CCS -j 1 --log-level FATAL
  $ samtools view -h $CRAMTMP/empty_fastq_xml_aligned_sorted_ccs.bam | grep -v "@PG"
  @HD\tVN:1.9\tSO:coordinate\tpb:3.0.7 (esc)
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)
