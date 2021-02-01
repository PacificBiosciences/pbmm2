  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align "$TESTDIR"/data/median.transcriptset.xml "$REF" "$CRAMTMP"/out_json_trans.json
  $ grep fileTypeId "$CRAMTMP"/out_json_trans.json | tr -d ' '
  "fileTypeId":"PacBio.DataSet.TranscriptAlignmentSet",
  $ grep path "$CRAMTMP"/out_json_trans.json | grep out_json_trans.transcriptalignmentset.xml | wc -l | tr -d ' '
  1
