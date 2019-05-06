  $ IN=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsorted.bam
  $ samtools view -H $CRAMTMP/unsorted.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsorted.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsorted.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsorted_pbi.bam --pbi
  $ samtools view -H $CRAMTMP/unsorted_pbi.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsorted_pbi.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsorted_pbi.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsorted_pbi.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsorted_pbi.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted.bam --sort
  $ samtools view -H $CRAMTMP/sorted.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sorted.bam.bai 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sorted.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted_nobai.bam --sort --no-bai
  $ samtools view -H $CRAMTMP/sorted_nobai.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sorted_nobai.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sorted_nobai.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sorted_nobai.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sorted_nobai.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted_pbi.bam --sort --pbi
  $ samtools view -H $CRAMTMP/sorted_pbi.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sorted_pbi.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sorted_pbi.bam.bai 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sorted_pbi.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/sorted_pbi.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsortedds.alignmentset.xml 2> $CRAMTMP/unsortedds.err || echo $?
  $ cut -f 8 -d '|' < $CRAMTMP/unsortedds.err
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H $CRAMTMP/unsortedds.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedds.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedds.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
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
  $ ls -alh $CRAMTMP/sortedds.bam.bai 2> /dev/null | wc -l | tr -d ' '
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
  $ ls -alh $CRAMTMP/unsortedjs.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh $CRAMTMP/unsortedjs.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/unsortedjs.json 2> /dev/null | wc -l | tr -d ' '
  1

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sortedjs.json --sort
  $ samtools view -H $CRAMTMP/sortedjs.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sortedjs.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh $CRAMTMP/sortedjs.bam.bai 2> /dev/null | wc -l | tr -d ' '
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
  $ ls -alh $CRAMTMP/unsortedoutstream.bam.bai 2> /dev/null | wc -l | tr -d ' '
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
  $ ls -alh $CRAMTMP/sortedoutstream.bam.bai 2> /dev/null | wc -l | tr -d ' '
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
  $ ls -alh $CRAMTMP/out_sub.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
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
  $ ls -alh $CRAMTMP/out_cons.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
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
  $ ls -alh $CRAMTMP/out_trans.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
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

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted.bam --sort -j 2 -J 2 -m 100M --log-level INFO 2>&1| grep INFO
  *Using 2 threads for alignments, 2 threads for sorting, and 200M bytes RAM for sorting. (glob)
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Merged sorted output from 0 files and 1 in-memory blocks (glob)
  *Generating BAI (glob)
  *Mapped Reads: 52 (glob)
  *Alignments: 96 (glob)
  *Mapped Bases: 242437 (glob)
  *Mean Mapped Concordance: 91* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Sort Merge Time: * (glob)
  *BAI Generation Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -j 500 2>&1| grep WARN
  *Requested more threads for alignment (500) than system-wide available* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail.bam -J 500 --sort 2>&1| grep AlignSettings | grep Requested
  *Requested more threads for sorting (500) and alignment (1) than system-wide available* (glob)
  *Requested more threads for sorting (500) than system-wide available* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/pass2.bam -j 1 -J 500 -m 500G
  *Requested 500 threads for sorting, without specifying --sort. Please check your input. (glob)
  *Requested 500G memory for sorting, without specifying --sort. Please check your input. (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail3.bam -j 1 -J 500 --sort --sort 2>&1| grep Requested
  *Requested more threads for sorting* (glob)
  *Requested more threads for sorting* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/fail4.bam -j 1 -J 2 --sort -m 1000G 2>&1
  *Trying to allocate more memory for sorting* (glob)
  [1]

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sort_percentage_4.bam -j 8 --sort --log-level INFO 2>&1 | grep "threads for alignments"
  *Using 6 threads for alignments, 2 threads for sorting, and 1.5G bytes RAM for sorting.* (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/default_parameters.bam --log-level DEBUG 2>&1| grep DEBUG
  *Minimap2 parameters* (glob)
  *Kmer size              : 19 (glob)
  *Minimizer window size  : 10 (glob)
  *Homopolymer compressed : true (glob)
  *Gap open 1             : 5 (glob)
  *Gap open 2             : 56 (glob)
  *Gap extension 1        : 4 (glob)
  *Gap extension 2        : 1 (glob)
  *Match score            : 2 (glob)
  *Mismatch penalty       : 5 (glob)
  *Z-drop                 : 400 (glob)
  *Z-drop inv             : 50 (glob)
  *Bandwidth              : 2000 (glob)
  *Long join flank ratio  : 0.5 (glob)

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/default_overrides.bam --log-level DEBUG -o 5 -O 56 -e 4 -E 1 -k 19 -w 10 -A 2 -B 5 -z 400 -Z 50 -r 2000 -L 0.4 2>&1| grep DEBUG
  *Minimap2 parameters* (glob)
  *Kmer size              : 19 (glob)
  *Minimizer window size  : 10 (glob)
  *Homopolymer compressed : true (glob)
  *Gap open 1             : 5 (glob)
  *Gap open 2             : 56 (glob)
  *Gap extension 1        : 4 (glob)
  *Gap extension 2        : 1 (glob)
  *Match score            : 2 (glob)
  *Mismatch penalty       : 5 (glob)
  *Z-drop                 : 400 (glob)
  *Z-drop inv             : 50 (glob)
  *Bandwidth              : 2000 (glob)
  *Long join flank ratio  : 0.4 (glob)

Test bam_sort
  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted_small.bam --sort -J 1 -m 1M --log-level TRACE --log-file $CRAMTMP/sorted_small.txt

Test that median filter does not fail
  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/median_output.bam --median-filter

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/bestn1.bam --best-n 1
  $ samtools view $CRAMTMP/bestn1.bam | wc -l | tr -d ' '
  52
