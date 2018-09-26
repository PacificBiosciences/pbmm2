  $ cd $TESTDIR/data
  $ IN=median.json
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ji_unsorted.bam
  $ samtools view -H $CRAMTMP/ji_unsorted.bam | grep "@HD" | grep "unknown" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_unsorted.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/ji_unsorted.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/ji_unsorted.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ji_sorted.bam --sort
  $ samtools view -H $CRAMTMP/ji_sorted.bam | grep "@HD" | grep "coordinate" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_sorted.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/ji_sorted.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/ji_sorted.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ji_unsortedds.alignmentset.xml 2> $CRAMTMP/ji_unsortedds.err || echo $?
  $ cut -f 8 -d '|' < $CRAMTMP/ji_unsortedds.err
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H $CRAMTMP/ji_unsortedds.bam | grep "@HD" | grep "unknown" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_unsortedds.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_unsortedds.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_unsortedds.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ji_sortedds.alignmentset.xml --sort 2> $CRAMTMP/ji_sortedds.err || echo $?
  $ cut -f 8 -d '|' < $CRAMTMP/ji_sortedds.err
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H $CRAMTMP/ji_sortedds.bam | grep "@HD" | grep "coordinate" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_sortedds.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_sortedds.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_sortedds.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ji_unsortedjs.json
  $ samtools view -H $CRAMTMP/ji_unsortedjs.bam | grep "@HD" | grep "unknown" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_unsortedjs.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_unsortedjs.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_unsortedjs.json 2> /dev/null | wc -l | sed 's/ //g'
  1

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/ji_sortedjs.json --sort
  $ samtools view -H $CRAMTMP/ji_sortedjs.bam | grep "@HD" | grep "coordinate" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_sortedjs.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_sortedjs.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_sortedjs.json 2> /dev/null | wc -l | sed 's/ //g'
  1

  $ $__PBTEST_PBMM2_EXE align $IN $REF > $CRAMTMP/ji_unsortedoutstream.bam
  $ samtools view -H $CRAMTMP/ji_unsortedoutstream.bam | grep "@HD" | grep "unknown" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_unsortedoutstream.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/ji_unsortedoutstream.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/ji_unsortedoutstream.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF > $CRAMTMP/ji_sortedoutstream.bam --sort
  $ samtools view -H $CRAMTMP/ji_sortedoutstream.bam | grep "@HD" | grep "coordinate" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/ji_sortedoutstream.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/ji_sortedoutstream.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/ji_sortedoutstream.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align median.json $REF $CRAMTMP/median_json_preset.bam --log-level INFO 2>&1 | grep "Setting to"
  *Setting to SUBREAD preset (glob)

  $ $__PBTEST_PBMM2_EXE align median.subreadset.json $REF $CRAMTMP/median_json_preset.bam --log-level INFO 2>&1 | grep "Setting to"
  *Setting to SUBREAD preset (glob)

  $ $__PBTEST_PBMM2_EXE align median.consensusreadset.json $REF $CRAMTMP/median_json_preset.bam --log-level INFO 2>&1 | grep "Setting to"
  *Setting to CCS preset (glob)

  $ $__PBTEST_PBMM2_EXE align median.transcriptset.json $REF $CRAMTMP/median_json_preset.bam --log-level INFO 2>&1 | grep "Setting to"
  *Setting to ISOSEQ preset (glob)
