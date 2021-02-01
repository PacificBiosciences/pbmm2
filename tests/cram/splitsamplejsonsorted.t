  $ cd "$TESTDIR"/data
  $ MERGED=merged.json
  $ REF=ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align "$MERGED" "$REF" "$CRAMTMP"/split_dataset_sorted_json.alignmentset.xml --split-by-sample --sort

  $ [[ -f "$CRAMTMP"/split_dataset_sorted_json.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.alignmentset.xml ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset_sorted_json.test_test.alignmentset.xml ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset_sorted_json.UCLA_1023.alignmentset.xml ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset_sorted_json.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset_sorted_json.test_test.bam ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset_sorted_json.UCLA_1023.bam ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset_sorted_json.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam.pbi ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset_sorted_json.test_test.bam.pbi ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset_sorted_json.UCLA_1023.bam.pbi ]] || echo "File does not exist!"

  $ grep "PacBio.AlignmentFile.AlignmentBamFile" "$CRAMTMP"/split_dataset_sorted_json.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.alignmentset.xml | tr -d '\t' | cut -f 4 -d ' ' | cut -f 2 -d '=' | tr -d '"' | sed 's|^.*/||g'
  split_dataset_sorted_json.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam
  $ grep "PacBio.AlignmentFile.AlignmentBamFile" "$CRAMTMP"/split_dataset_sorted_json.test_test.alignmentset.xml | tr -d '\t' | cut -f 4 -d ' ' | cut -f 2 -d '=' | tr -d '"' | sed 's|^.*/||g'
  split_dataset_sorted_json.test_test.bam
  $ grep "PacBio.AlignmentFile.AlignmentBamFile" "$CRAMTMP"/split_dataset_sorted_json.UCLA_1023.alignmentset.xml | tr -d '\t' | cut -f 4 -d ' ' | cut -f 2 -d '=' | tr -d '"' | sed 's|^.*/||g'
  split_dataset_sorted_json.UCLA_1023.bam

  $ grep -c split_dataset_sorted_json.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.alignmentset.xml "$CRAMTMP"/split_dataset_sorted_json.json
  1
  $ grep -c split_dataset_sorted_json.test_test.alignmentset.xml "$CRAMTMP"/split_dataset_sorted_json.json
  1
  $ grep -c split_dataset_sorted_json.UCLA_1023.alignmentset.xml "$CRAMTMP"/split_dataset_sorted_json.json
  1

  $ ID=$("$SAMTOOLS" view -F 4 -H "$CRAMTMP"/split_dataset_sorted_json.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam | grep 3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494 | grep -v "@PG	ID:samtools" | cut -f 2 | cut -f 2 -d ':')
  $ "$SAMTOOLS" view -F 4 "$CRAMTMP"/split_dataset_sorted_json.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam | grep -vc ${ID} | grep -v "@PG	ID:samtools" | tr -d ' '
  0
  $ "$SAMTOOLS" view -F 4 "$CRAMTMP"/split_dataset_sorted_json.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam | grep -v "@PG	ID:samtools" | wc -l | tr -d ' '
  9

  $ ID=$("$SAMTOOLS" view -F 4 -H "$CRAMTMP"/split_dataset_sorted_json.test_test.bam | grep test_test | grep -v "@PG	ID:samtools" | cut -f 2 | cut -f 2 -d ':')
  $ "$SAMTOOLS" view -F 4 "$CRAMTMP"/split_dataset_sorted_json.test_test.bam | grep -vc ${ID} | grep -v "@PG	ID:samtools" | tr -d ' '
  0
  $ "$SAMTOOLS" view -F 4 "$CRAMTMP"/split_dataset_sorted_json.test_test.bam | grep -v "@PG	ID:samtools" | wc -l | tr -d ' '
  11

  $ ID=$("$SAMTOOLS" view -F 4 -H "$CRAMTMP"/split_dataset_sorted_json.UCLA_1023.bam | grep UCLA_1023 | grep -v "@PG	ID:samtools" | cut -f 2 | cut -f 2 -d ':')
  $ "$SAMTOOLS" view -F 4 "$CRAMTMP"/split_dataset_sorted_json.UCLA_1023.bam | grep -vc ${ID} | grep -v "@PG	ID:samtools" | tr -d ' '
  0
  $ "$SAMTOOLS" view -F 4 "$CRAMTMP"/split_dataset_sorted_json.UCLA_1023.bam | grep -v "@PG	ID:samtools" | wc -l | tr -d ' '
  10
