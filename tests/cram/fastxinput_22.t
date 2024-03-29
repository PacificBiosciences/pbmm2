  $ BAM="$TESTDIR"/data/m54075_180905_225130.ccs.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ samtools view "$BAM" | awk '{ print ">"$1"\n"$10 }' > "$CRAMTMP"/m54075_180905_225130.fasta
  $ FASTA="$CRAMTMP"/m54075_180905_225130.fasta

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_unsorted.bam --preset CCS 2>&1
  *Input is FASTA.* (glob)
  $ samtools view -H "$CRAMTMP"/fasta_unsorted.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsorted.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsorted.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsorted.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_sorted.bam --preset CCS --sort
  *Input is FASTA.* (glob)
  $ samtools view -H "$CRAMTMP"/fasta_sorted.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_sorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_sorted.bam.bai 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_sorted.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_sorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_sorted.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_unsortedds.alignmentset.xml --preset CCS 2> "$CRAMTMP"/fasta_unsortedds.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/fasta_unsortedds.err
  *Input is FASTA.* (glob)
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/fasta_unsortedds.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedds.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedds.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedds.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedds.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedds.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_sortedds.alignmentset.xml --preset CCS --sort 2> "$CRAMTMP"/fasta_sortedds.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/fasta_sortedds.err
  *Input is FASTA.* (glob)
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/fasta_sortedds.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_sortedds.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_sortedds.bam.bai 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_sortedds.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_sortedds.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_sortedds.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_unsortedjs.json --preset CCS
  *Input is FASTA.* (glob)
  *Unsupported input type* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_sortedjs.json --preset CCS --sort
  *Input is FASTA.* (glob)
  *Unsupported input type* (glob)
  [1]

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" > "$CRAMTMP"/fasta_unsortedoutstream.bam --preset CCS
  *Input is FASTA.* (glob)
  $ samtools view -H "$CRAMTMP"/fasta_unsortedoutstream.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedoutstream.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" > "$CRAMTMP"/fasta_sortedoutstream.bam --preset CCS --sort
  *Input is FASTA.* (glob)
  $ samtools view -H "$CRAMTMP"/fasta_sortedoutstream.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_sortedoutstream.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_sortedoutstream.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_sortedoutstream.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_sortedoutstream.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_sortedoutstream.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_unsortedccs.consensusalignmentset.xml --preset CCS 2> "$CRAMTMP"/fasta_unsortedccs.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/fasta_unsortedccs.err
  *Input is FASTA.* (glob)
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/fasta_unsortedccs.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedccs.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align -j 1 "$REF" "$FASTA" "$CRAMTMP"/fasta_unsortedts.transcriptalignmentset.xml --preset CCS 2> "$CRAMTMP"/fasta_unsortedts.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/fasta_unsortedts.err
  *Input is FASTA.* (glob)
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/fasta_unsortedts.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedts.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedts.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedts.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/fasta_unsortedts.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/fasta_unsortedts.json 2> /dev/null | wc -l | tr -d ' '
  0

  $ "$PBMM2" align "$REF" "$FASTA" "$CRAMTMP"/sorted_fasta_verbose.bam --preset CCS --sort -j 2 -J 2 -m 100M --log-level INFO 2>&1
  *Using 2 threads for alignments, 2 threads for sorting, and 200M bytes RAM for sorting. (glob)
  *Input is FASTA.* (glob)
  *READ input file: *m54075_180905_225130.fasta* (glob)
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
