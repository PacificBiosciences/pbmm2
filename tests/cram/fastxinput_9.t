  $ BAM="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ samtools view "$BAM" | awk '{ print ">"$1"\n"$10 }' > "$CRAMTMP"/median.fasta
  $ FASTA="$CRAMTMP"/median.fasta

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_unsortedccs.consensusalignmentset.xml --preset SUBREAD 2> "$CRAMTMP"/fasta_unsortedccs.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/fasta_unsortedccs.err
  *Input is FASTA.* (glob)
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/fasta_unsortedccs.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.json 2> /dev/null | wc -l | tr -d ' '
  0
