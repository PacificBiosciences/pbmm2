  $ IN=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsorted.bam
  $ samtools view -H $CRAMTMP/unsorted.bam | grep "@HD" | grep "unknown" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/unsorted.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/unsorted.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/unsorted.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted.bam --sort
  $ samtools view -H $CRAMTMP/sorted.bam | grep "@HD" | grep "coordinate" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/sorted.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/sorted.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/sorted.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsortedds.alignmentset.xml
  $ samtools view -H $CRAMTMP/unsortedds.bam | grep "@HD" | grep "unknown" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/unsortedds.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/unsortedds.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/unsortedds.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sortedds.alignmentset.xml --sort
  $ samtools view -H $CRAMTMP/sortedds.bam | grep "@HD" | grep "coordinate" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/sortedds.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/sortedds.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/sortedds.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsortedjs.json
  $ samtools view -H $CRAMTMP/unsortedjs.bam | grep "@HD" | grep "unknown" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/unsortedjs.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/unsortedjs.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/unsortedjs.json 2> /dev/null | wc -l | sed 's/ //g'
  1

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sortedjs.json --sort
  $ samtools view -H $CRAMTMP/sortedjs.bam | grep "@HD" | grep "coordinate" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/sortedjs.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/sortedjs.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/sortedjs.json 2> /dev/null | wc -l | sed 's/ //g'
  1

  $ $__PBTEST_PBMM2_EXE align $IN $REF > $CRAMTMP/unsortedoutstream.bam
  $ samtools view -H $CRAMTMP/unsortedoutstream.bam | grep "@HD" | grep "unknown" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/unsortedoutstream.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/unsortedoutstream.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/unsortedoutstream.json 2> /dev/null | wc -l | sed 's/ //g'
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF > $CRAMTMP/sortedoutstream.bam --sort
  $ samtools view -H $CRAMTMP/sortedoutstream.bam | grep "@HD" | grep "coordinate" | wc -l | sed 's/ //g'
  1
  $ ls -alh $CRAMTMP/sortedoutstream.bam.pbi 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/sortedoutstream.*.xml 2> /dev/null | wc -l | sed 's/ //g'
  0
  $ ls -alh $CRAMTMP/sortedoutstream.json 2> /dev/null | wc -l | sed 's/ //g'
  0

