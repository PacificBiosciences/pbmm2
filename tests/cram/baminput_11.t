  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/unsortedjs.json
  $ "$SAMTOOLS" view -H "$CRAMTMP"/unsortedjs.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/unsortedjs.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/unsortedjs.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsortedjs.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsortedjs.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/unsortedjs.json 2> /dev/null | wc -l | tr -d ' '
  1
