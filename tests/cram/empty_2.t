  $ FASTQ="$TESTDIR"/data/empty.fastq
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta
  $ REFXML="$TESTDIR"/data/ecoli.referenceset.xml

FASTQ tests
  $ "$PBMM2" align "$REF" "$FASTQ" "$CRAMTMP"/empty_fastq_aligned_unsorted.bam -j 1 --log-level FATAL
  $ "$SAMTOOLS" view -h "$CRAMTMP"/empty_fastq_aligned_unsorted.bam | grep -v "@PG" | grep -v "@HD"
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ "$PBMM2" align "$REF" "$FASTQ" "$CRAMTMP"/empty_fastq_aligned_sorted.bam --sort -j 1 --log-level FATAL
  $ "$SAMTOOLS" view -h "$CRAMTMP"/empty_fastq_aligned_sorted.bam | grep -v "@PG" | grep -v "@HD"
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ "$PBMM2" align "$REF" "$FASTQ" "$CRAMTMP"/empty_fastq_aligned_sorted_ccs.bam --sort --preset CCS -j 1 --log-level FATAL
  $ "$SAMTOOLS" view -h "$CRAMTMP"/empty_fastq_aligned_sorted_ccs.bam | grep -v "@PG" | grep -v "@HD"
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ "$PBMM2" align "$REFXML" "$FASTQ" "$CRAMTMP"/empty_fastq_xml_aligned_unsorted.bam -j 1 --log-level FATAL
  $ "$SAMTOOLS" view -h "$CRAMTMP"/empty_fastq_xml_aligned_unsorted.bam | grep -v "@PG" | grep -v "@HD"
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ "$PBMM2" align "$REFXML" "$FASTQ" "$CRAMTMP"/empty_fastq_xml_aligned_sorted.bam --sort -j 1 --log-level FATAL
  $ "$SAMTOOLS" view -h "$CRAMTMP"/empty_fastq_xml_aligned_sorted.bam | grep -v "@PG" | grep -v "@HD"
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)

  $ "$PBMM2" align "$REFXML" "$FASTQ" "$CRAMTMP"/empty_fastq_xml_aligned_sorted_ccs.bam --sort --preset CCS -j 1 --log-level FATAL
  $ "$SAMTOOLS" view -h "$CRAMTMP"/empty_fastq_xml_aligned_sorted_ccs.bam | grep -v "@PG" | grep -v "@HD"
  @SQ\tSN:ecoliK12_pbi_March2013\tLN:4642522 (esc)
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)
