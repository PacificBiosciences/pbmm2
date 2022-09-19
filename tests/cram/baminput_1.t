  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/unsorted.bam
  $ samtools view -H "$CRAMTMP"/unsorted.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/unsorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsorted.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsorted.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsorted.json 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsorted.ref.collapsed.fasta 2> /dev/null | wc -l | tr -d ' '
  0
