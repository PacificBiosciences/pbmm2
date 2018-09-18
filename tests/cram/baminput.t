  $ IN=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsorted.bam
  $ samtools view -H $CRAMTMP/unsorted.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsorted.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted.bam --sort
  $ samtools view -H $CRAMTMP/sorted.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sorted.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsortedds.alignmentset.xml 2> $CRAMTMP/unsortedds.err || echo $?
  $ cut -f 8 -d '|' < $CRAMTMP/unsortedds.err
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H $CRAMTMP/unsortedds.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedds.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedds.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedds.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sortedds.alignmentset.xml --sort 2> $CRAMTMP/sortedds.err || echo $?
  $ cut -f 8 -d '|' < $CRAMTMP/sortedds.err
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H $CRAMTMP/sortedds.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sortedds.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sortedds.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sortedds.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsortedjs.json
  $ samtools view -H $CRAMTMP/unsortedjs.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedjs.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedjs.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedjs.json 2> /dev/null | wc -l | tr -d ' '
  1

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sortedjs.json --sort
  $ samtools view -H $CRAMTMP/sortedjs.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sortedjs.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sortedjs.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sortedjs.json 2> /dev/null | wc -l | tr -d ' '
  1

  $ $__PBTEST_PBMM2_EXE align $IN $REF > $CRAMTMP/unsortedoutstream.bam
  $ samtools view -H $CRAMTMP/unsortedoutstream.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedoutstream.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsortedoutstream.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsortedoutstream.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF > $CRAMTMP/sortedoutstream.bam --sort
  $ samtools view -H $CRAMTMP/sortedoutstream.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sortedoutstream.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sortedoutstream.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sortedoutstream.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.subreadset.xml $REF $CRAMTMP/out_sub.subreadset.xml 2> $CRAMTMP/out_sub.err || echo $?
  1
  $ cut -f 8 -d '|' < $CRAMTMP/out_sub.err
  - Output has to be an alignment dataset! Please use alignmentset.xml, consensusalignmentset.xml, or transcriptalignmentset.xml!

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.subreadset.xml $REF $CRAMTMP/out_sub.alignmentset.xml
  $ ls -alh $CRAMTMP/out_sub.bam 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/out_sub.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/out_sub.alignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/out_sub.consensusalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/out_sub.transcriptalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/out_sub.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.consensusreadset.xml $REF $CRAMTMP/out_cons.consensusalignmentset.xml
  $ ls -alh $CRAMTMP/out_cons.bam 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/out_cons.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/out_cons.alignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/out_cons.consensusalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/out_cons.transcriptalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/out_cons.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.transcriptset.xml $REF $CRAMTMP/out_trans.transcriptalignmentset.xml
  $ echo $?
  0
  $ ls -alh $CRAMTMP/out_trans.bam 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/out_trans.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/out_trans.alignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/out_trans.consensusalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/out_trans.transcriptalignmentset.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/out_trans.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.subreadset.xml $REF $CRAMTMP/out_cons_fail.consensusalignmentset.xml 2> $CRAMTMP/out_cons_fail.err || echo $?
  1
  $ cut -f 8 -d '|' < $CRAMTMP/out_cons_fail.err
  - Unsupported dataset combination! Input SubreadSet with output ConsensusReadSet! Please use AlignmentSet as output XML type!

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.consensusreadset.xml $REF $CRAMTMP/out_trans_fail.transcriptalignmentset.xml 2> $CRAMTMP/out_trans_fail.err || echo $?
  1
  $ cut -f 8 -d '|' < $CRAMTMP/out_trans_fail.err
  - Unsupported dataset combination! Input ConsensusReadSet with output TranscriptSet! Please use ConsensusAlignmentSet as output XML type!

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.transcriptset.xml $REF $CRAMTMP/out_align_fail.alignmentset.xml 2> $CRAMTMP/out_align_fail.err || echo $?
  1
  $ cut -f 8 -d '|' < $CRAMTMP/out_align_fail.err
  - Unsupported dataset combination! Input TranscriptSet with output AlignmentSet! Please use TranscriptAlignmentSet as output XML type!

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.subreadset.xml $REF $CRAMTMP/out_json_sub.json
  $ grep fileTypeId $CRAMTMP/out_json_sub.json | tr -d ' '
  "fileTypeId":"PacBio.DataSet.AlignmentSet",
  $ grep path $CRAMTMP/out_json_sub.json | grep out_json_sub.alignmentset.xml | wc -l | tr -d ' '
  1

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.consensusreadset.xml $REF $CRAMTMP/out_json_ccs.json
  $ grep fileTypeId $CRAMTMP/out_json_ccs.json | tr -d ' '
  "fileTypeId":"PacBio.DataSet.ConsensusAlignmentSet",
  $ grep path $CRAMTMP/out_json_ccs.json | grep out_json_ccs.consensusalignmentset.xml | wc -l | tr -d ' '
  1

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.transcriptset.xml $REF $CRAMTMP/out_json_trans.json
  $ grep fileTypeId $CRAMTMP/out_json_trans.json | tr -d ' '
  "fileTypeId":"PacBio.DataSet.TranscriptAlignmentSet",
  $ grep path $CRAMTMP/out_json_trans.json | grep out_json_trans.transcriptalignmentset.xml | wc -l | tr -d ' '
  1

  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.transcriptset.xml $REF $CRAMTMP/out_trans_upper.TranscriptAlignmentSet.XmL || echo $?
  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.transcriptset.xml $REF $CRAMTMP/out_json_upper.JsON || echo $?
  $ $__PBTEST_PBMM2_EXE align $TESTDIR/data/median.transcriptset.xml $REF $CRAMTMP/out_xml_upper.XML 2> $CRAMTMP/out_xml_upper.err  || echo $?
  1
  $ cut -f 8 -d '|' < $CRAMTMP/out_xml_upper.err
  - Output is XML, but of unknown type! Please use alignmentset.xml, consensusalignmentset.xml, or transcriptalignmentset.xml
