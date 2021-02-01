  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align "$TESTDIR"/data/median.subreadset.xml "$REF" "$CRAMTMP"/out_cons_fail.consensusalignmentset.xml 2> "$CRAMTMP"/out_cons_fail.err || echo $?
  1
  $ cut -f 8 -d '|' < "$CRAMTMP"/out_cons_fail.err
  *Unsupported dataset combination! Input SubreadSet with output ConsensusReadSet! Please use AlignmentSet as output XML type!* (glob)

  $ "$PBMM2" align "$TESTDIR"/data/median.consensusreadset.xml "$REF" "$CRAMTMP"/out_trans_fail.transcriptalignmentset.xml 2> "$CRAMTMP"/out_trans_fail.err || echo $?
  1
  $ cut -f 8 -d '|' < "$CRAMTMP"/out_trans_fail.err
  *Unsupported dataset combination! Input ConsensusReadSet with output TranscriptSet! Please use ConsensusAlignmentSet as output XML type!* (glob)

  $ "$PBMM2" align "$TESTDIR"/data/median.transcriptset.xml "$REF" "$CRAMTMP"/out_align_fail.alignmentset.xml 2> "$CRAMTMP"/out_align_fail.err || echo $?
  1
  $ cut -f 8 -d '|' < "$CRAMTMP"/out_align_fail.err
  *Unsupported dataset combination! Input TranscriptSet with output AlignmentSet! Please use TranscriptAlignmentSet as output XML type!* (glob)

  $ "$PBMM2" align "$TESTDIR"/data/median.subreadset.xml "$REF" "$CRAMTMP"/out_json_sub.json
  $ grep fileTypeId "$CRAMTMP"/out_json_sub.json | tr -d ' '
  "fileTypeId":"PacBio.DataSet.AlignmentSet",
  $ grep path "$CRAMTMP"/out_json_sub.json | grep out_json_sub.alignmentset.xml | wc -l | tr -d ' '
  1
