  $ BAM="$TESTDIR"/data/m54075_180905_225130.ccs.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ samtools view "$BAM" | awk '{ print "@"$1"\n"$10"\n+\n"$11 }' > "$CRAMTMP"/m54075_180905_225130.fastq
  $ FASTQ="$CRAMTMP"/m54075_180905_225130.fastq

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/fastq_unsorted.bam --preset CCS 2>&1
  *Input is FASTQ.* (glob)
  $ samtools view -H "$CRAMTMP"/fastq_unsorted.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsorted.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsorted.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsorted.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/fastq_sorted.bam --preset CCS --sort
  *Input is FASTQ.* (glob)
  $ samtools view -H "$CRAMTMP"/fastq_sorted.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_sorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sorted.bam.bai 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_sorted.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sorted.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/fastq_unsortedds.alignmentset.xml --preset CCS 2> "$CRAMTMP"/fastq_unsortedds.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/fastq_unsortedds.err
  *Input is FASTQ.* (glob)
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/fastq_unsortedds.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedds.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedds.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsortedds.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsortedds.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedds.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/fastq_sortedds.alignmentset.xml --preset CCS --sort 2> "$CRAMTMP"/fastq_sortedds.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/fastq_sortedds.err
  *Input is FASTQ.* (glob)
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/fastq_sortedds.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_sortedds.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_sortedds.bam.bai 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_sortedds.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sortedds.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_sortedds.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/fastq_unsortedjs.json --preset CCS
  *Input is FASTQ.* (glob)
  *Unsupported input type* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/fastq_sortedjs.json --preset CCS --sort
  *Input is FASTQ.* (glob)
  *Unsupported input type* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" > "$CRAMTMP"/fastq_unsortedoutstream.bam --preset CCS
  *Input is FASTQ.* (glob)
  $ samtools view -H "$CRAMTMP"/fastq_unsortedoutstream.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedoutstream.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsortedoutstream.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsortedoutstream.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsortedoutstream.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsortedoutstream.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" > "$CRAMTMP"/fastq_sortedoutstream.bam --preset CCS --sort
  *Input is FASTQ.* (glob)
  $ samtools view -H "$CRAMTMP"/fastq_sortedoutstream.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_sortedoutstream.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/fastq_unsortedccs.consensusalignmentset.xml --preset CCS 2> "$CRAMTMP"/fastq_unsortedccs.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/fastq_unsortedccs.err
  *Input is FASTQ.* (glob)
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/fastq_unsortedccs.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedccs.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedccs.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsortedccs.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fastq_unsortedccs.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedccs.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTQ" "$CRAMTMP"/fastq_unsortedts.transcriptalignmentset.xml --preset CCS 2> "$CRAMTMP"/fastq_unsortedts.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/fastq_unsortedts.err
  *Input is FASTQ.* (glob)
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/fastq_unsortedts.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedts.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedts.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedts.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fastq_unsortedts.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align "$REF" "$FASTQ" "$CRAMTMP"/sorted_fastq_verbose.bam --preset CCS --sort -j 2 -J 2 -m 100M --log-level INFO 2>&1
  *Using 2 threads for alignments, 2 threads for sorting, and 200M bytes RAM for sorting. (glob)
  *Input is FASTQ.* (glob)
  *READ input file: *m54075_180905_225130.fastq* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Merged sorted output from 0 files and 1 in-memory blocks (glob)
  *Generating BAI (glob)
  *Mapped Reads: 7 (glob)
  *Alignments: 7 (glob)
  *Mapped Bases: 11994 (glob)
  *Mean Gap-Compressed Sequence Identity* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Sort Merge Time: * (glob)
  *BAI Generation Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)
