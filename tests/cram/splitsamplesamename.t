  $ MERGED=$TESTDIR/data/merged.same.dataset.xml
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $MERGED $REF $CRAMTMP/splitname.bam --split-by-sample

  $ ls -l $CRAMTMP/splitname.*.bam | wc -l | tr -d ' '
  3

  $ ID=$(samtools view -F 4 -H $CRAMTMP/splitname.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam | grep 3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494 | cut -f 2 | cut -f 2 -d ':')
  $ samtools view -F 4 $CRAMTMP/splitname.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam | grep -vc ${ID} | tr -d ' '
  0
  $ samtools view -F 4 $CRAMTMP/splitname.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam | wc -l | tr -d ' '
  9

  $ ID=$(samtools view -F 4 -H $CRAMTMP/splitname.test-0.bam | grep "test(" | cut -f 2 | cut -f 2 -d ':')
  $ samtools view -F 4 $CRAMTMP/splitname.test-0.bam | grep -vc ${ID} | tr -d ' '
  0
  $ samtools view -F 4 $CRAMTMP/splitname.test-0.bam | wc -l | tr -d ' '
  10

  $ ID=$(samtools view -F 4 -H $CRAMTMP/splitname.test-1.bam | grep "test)" | cut -f 2 | cut -f 2 -d ':')
  $ samtools view -F 4 $CRAMTMP/splitname.test-1.bam | grep -vc ${ID} | tr -d ' '
  0
  $ samtools view -F 4 $CRAMTMP/splitname.test-1.bam | wc -l | tr -d ' '
  11

  $ $__PBTEST_PBMM2_EXE align $MERGED $REF $CRAMTMP/split_dataset.alignmentset.xml --split-by-sample

  $ [[ -f $CRAMTMP/split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.alignmentset.xml ]] || echo "File does not exist!"
  $ [[ -f $CRAMTMP/split_dataset.test-0.alignmentset.xml ]] || echo "File does not exist!"
  $ [[ -f $CRAMTMP/split_dataset.test-1.alignmentset.xml ]] || echo "File does not exist!"
  $ [[ -f $CRAMTMP/split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam ]] || echo "File does not exist!"
  $ [[ -f $CRAMTMP/split_dataset.test-0.bam ]] || echo "File does not exist!"
  $ [[ -f $CRAMTMP/split_dataset.test-1.bam ]] || echo "File does not exist!"
  $ [[ -f $CRAMTMP/split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam.pbi ]] || echo "File does not exist!"
  $ [[ -f $CRAMTMP/split_dataset.test-0.bam.pbi ]] || echo "File does not exist!"
  $ [[ -f $CRAMTMP/split_dataset.test-1.bam.pbi ]] || echo "File does not exist!"

  $ grep "PacBio.AlignmentFile.AlignmentBamFile" $CRAMTMP/split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.alignmentset.xml | tr -d '\t' | cut -f 3 -d ' ' | cut -f 2 -d '=' | tr -d '"'
  split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.bam
  $ grep "PacBio.AlignmentFile.AlignmentBamFile" $CRAMTMP/split_dataset.test-0.alignmentset.xml | tr -d '\t' | cut -f 3 -d ' ' | cut -f 2 -d '=' | tr -d '"'
  split_dataset.test-0.bam
  $ grep "PacBio.AlignmentFile.AlignmentBamFile" $CRAMTMP/split_dataset.test-1.alignmentset.xml | tr -d '\t' | cut -f 3 -d ' ' | cut -f 2 -d '=' | tr -d '"'
  split_dataset.test-1.bam

  $ grep -c split_dataset.3260208_188nM-GTAC_2xGCratio_LP7_100fps_15min_5kEColi_SP2p1_3uMSSB_BA243494.alignmentset.xml $CRAMTMP/split_dataset.json
  1
  $ grep -c split_dataset.test-0.alignmentset.xml $CRAMTMP/split_dataset.json
  1
  $ grep -c split_dataset.test-1.alignmentset.xml $CRAMTMP/split_dataset.json
  1
