  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align "$TESTDIR"/data/median.subreadset.xml "$REF" "$CRAMTMP"/out_sub.subreadset.xml 2> "$CRAMTMP"/out_sub.err || echo $?
  1
  $ cut -f 8 -d '|' < "$CRAMTMP"/out_sub.err
  *Output has to be an alignment dataset! Please use alignmentset.xml, consensusalignmentset.xml, or transcriptalignmentset.xml!* (glob)

  $ "$PBMM2" align "$TESTDIR"/data/median.subreadset.xml "$REF" "$CRAMTMP"/out_sub.alignmentset.xml
  $ ls -alh "$CRAMTMP"/out_sub.bam 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/out_sub.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/out_sub.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_sub.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_sub.alignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/out_sub.consensusalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_sub.transcriptalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_sub.json 2> /dev/null | wc -l | tr -d ' '
  0
