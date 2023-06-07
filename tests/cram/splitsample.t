  $ MERGED="$TESTDIR"/data/merged.dataset.xml
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta
  $ NO_SM_BIOSAMPLES="$TESTDIR"/data/no_sm_biosamples.subreadset.xml

  $ "$PBMM2" align -j 1 "$MERGED" "$REF" "$CRAMTMP"/split.bam --preset SUBREAD --split-by-sample
  $ ls -l "$CRAMTMP"/split.*.bam | wc -l | tr -d ' '
  3

  $ ID=$(samtools view -F 4 -H "$CRAMTMP"/split.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam | grep 3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494 | grep -v "@PG	ID:samtools" | cut -f 2 | cut -f 2 -d ':')
  $ samtools view -F 4 "$CRAMTMP"/split.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam | grep -vc ${ID} | grep -v "@PG	ID:samtools" | tr -d ' '
  0
  $ samtools view -F 4 "$CRAMTMP"/split.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam | grep -v "@PG	ID:samtools" | wc -l | tr -d ' '
  9

  $ ID=$(samtools view -F 4 -H "$CRAMTMP"/split.test_test.bam | grep "test test" | grep -v "@PG	ID:samtools" | cut -f 2 | cut -f 2 -d ':')
  $ samtools view -F 4 "$CRAMTMP"/split.test_test.bam | grep -vc ${ID} | grep -v "@PG	ID:samtools" | tr -d ' '
  0
  $ samtools view -F 4 "$CRAMTMP"/split.test_test.bam | grep -v "@PG	ID:samtools" | wc -l | tr -d ' '
  11

  $ ID=$(samtools view -F 4 -H "$CRAMTMP"/split.UCLA_1023.bam | grep "UCLA 1023" | grep -v "@PG	ID:samtools" | cut -f 2 | cut -f 2 -d ':')
  $ samtools view -F 4 "$CRAMTMP"/split.UCLA_1023.bam | grep -vc ${ID} | grep -v "@PG	ID:samtools" | tr -d ' '
  0
  $ samtools view -F 4 "$CRAMTMP"/split.UCLA_1023.bam | grep -v "@PG	ID:samtools" | wc -l | tr -d ' '
  10

  $ "$PBMM2" align -j 1 "$MERGED" "$REF" "$CRAMTMP"/split_dataset.alignmentset.xml --preset SUBREAD --split-by-sample

  $ [[ -f "$CRAMTMP"/split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.alignmentset.xml ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset.test_test.alignmentset.xml ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset.UCLA_1023.alignmentset.xml ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset.test_test.bam ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset.UCLA_1023.bam ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam.pbi ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset.test_test.bam.pbi ]] || echo "File does not exist!"
  $ [[ -f "$CRAMTMP"/split_dataset.UCLA_1023.bam.pbi ]] || echo "File does not exist!"

  $ grep "PacBio.AlignmentFile.AlignmentBamFile" "$CRAMTMP"/split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.alignmentset.xml | tr -d '\t' | cut -f 4 -d ' ' | cut -f 2 -d '=' | tr -d '"' | sed 's|^.*/||g'
  split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam
  $ grep "PacBio.AlignmentFile.AlignmentBamFile" "$CRAMTMP"/split_dataset.test_test.alignmentset.xml | tr -d '\t' | cut -f 4 -d ' ' | cut -f 2 -d '=' | tr -d '"' | sed 's|^.*/||g'
  split_dataset.test_test.bam
  $ grep "PacBio.AlignmentFile.AlignmentBamFile" "$CRAMTMP"/split_dataset.UCLA_1023.alignmentset.xml | tr -d '\t' | cut -f 4 -d ' ' | cut -f 2 -d '=' | tr -d '"' | sed 's|^.*/||g'
  split_dataset.UCLA_1023.bam

  $ grep -c split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.alignmentset.xml "$CRAMTMP"/split_dataset.json
  1
  $ grep -c split_dataset.test_test.alignmentset.xml "$CRAMTMP"/split_dataset.json
  1
  $ grep -c split_dataset.UCLA_1023.alignmentset.xml "$CRAMTMP"/split_dataset.json
  1

When both --split-by-sample and --sample were set, expect to see only one bam file with SM overridden.
  $ IN="$TESTDIR"/data/merged.consensusreadset.xml
  $ "$PBMM2" align -j 1 "$REF" "$IN" "$CRAMTMP"/splitsampleoverride.consensusalignmentset.xml --preset SUBREAD --sort -j 8 --split-by-sample --sample "MySample"  2>&1 | grep -F -v 'Requested more threads'
  *Options --split-by-sample and --sample are mutually exclusive. Option --sample will be applied and --split-by-sample is ignored! (glob)
  *Offending bio sample names. BAM contains 'bamSample' and XML contains 'UCLA 1023'. Will ignore XML bio sample name.* (glob)
  $ [[ -f "$CRAMTMP"/splitsampleoverride.bam ]] || echo "File does not exist!"
  $ samtools view -H "$CRAMTMP"/splitsampleoverride.bam | grep "@RG" | grep -v "@PG	ID:samtools" | cut -f 6 | sort | uniq
  SM:MySample

  $ "$PBMM2" align -j 1 --split-by-sample "$NO_SM_BIOSAMPLES" "$REF" "$CRAMTMP"/split-no-sm.bam --preset SUBREAD
  $ samtools view -H "$CRAMTMP"/split-no-sm.UnnamedSample.bam | grep -c "@RG"
  1
  $ samtools view -H "$CRAMTMP"/split-no-sm.UnnamedSample.bam | grep -c "SM:UnnamedSample"
  1
