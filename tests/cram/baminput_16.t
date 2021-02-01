  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$TESTDIR"/data/median.consensusreadset.xml "$REF" "$CRAMTMP"/out_cons.consensusalignmentset.xml
  $ ls -alh "$CRAMTMP"/out_cons.bam 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/out_cons.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/out_cons.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_cons.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_cons.alignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_cons.consensusalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/out_cons.transcriptalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_cons.json 2> /dev/null | wc -l | tr -d ' '
  0
