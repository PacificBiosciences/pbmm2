  $ BAM="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ samtools view "$BAM" | awk '{ print ">"$1"\n"$10 }' > "$CRAMTMP"/median.fasta
  $ FASTA="$CRAMTMP"/median.fasta

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_unsortedjs.json --preset SUBREAD
  *Input is FASTA.* (glob)
  *Unsupported input type* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_sortedjs.json --preset SUBREAD --sort
  *Input is FASTA.* (glob)
  *Unsupported input type* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" --preset SUBREAD > "$CRAMTMP"/fasta_unsortedoutstream.bam
  *Input is FASTA.* (glob)
  $ samtools view -H "$CRAMTMP"/fasta_unsortedoutstream.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.json 2> /dev/null | wc -l | tr -d ' '
  0
