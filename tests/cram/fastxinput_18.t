  $ BAM="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ samtools view "$BAM" | awk '{ print "@"$1"\n"$10"\n+\n"$11 }' > "$CRAMTMP"/median.fastq
  $ FASTQ="$CRAMTMP"/median.fastq

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" --preset SUBREAD > "$CRAMTMP"/fastq_sortedoutstream.bam --sort
  *Input is FASTQ.* (glob)
  $ samtools view -H "$CRAMTMP"/fastq_sortedoutstream.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.json 2> /dev/null | wc -l | tr -d ' '
  0
