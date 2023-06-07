  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$TESTDIR"/data/median.consensusreadset.xml "$REF" "$CRAMTMP"/out_json_ccs.json --preset CCS
  $ grep fileTypeId "$CRAMTMP"/out_json_ccs.json | tr -d ' '
  "fileTypeId":"PacBio.DataSet.ConsensusAlignmentSet",
  $ grep path "$CRAMTMP"/out_json_ccs.json | grep out_json_ccs.consensusalignmentset.xml | wc -l | tr -d ' '
  1
