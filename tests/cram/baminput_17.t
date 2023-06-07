  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$TESTDIR"/data/median.transcriptset.xml "$REF" "$CRAMTMP"/out_trans.transcriptalignmentset.xml --preset ISOSEQ
  $ echo $?
  0
  $ ls -alh "$CRAMTMP"/out_trans.bam 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/out_trans.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/out_trans.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_trans.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_trans.alignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_trans.consensusalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/out_trans.transcriptalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/out_trans.json 2> /dev/null | wc -l | tr -d ' '
  0
