  $ BAM="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ samtools view "$BAM" | awk '{ print ">"$1"\n"$10 }' > "$CRAMTMP"/median.fasta
  $ FASTA="$CRAMTMP"/median.fasta

  $ "$PBMM2" align "$REF" "$FASTA" "$CRAMTMP"/sorted_fasta_verbose.bam --preset SUBREAD --sort -j 2 -J 2 -m 100M --log-level INFO 2>&1
  *Using 2 threads for alignments, 2 threads for sorting, and 200M bytes RAM for sorting. (glob)
  *Input is FASTA.* (glob)
  *READ input file: *median.fasta* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Merged sorted output from 0 files and 1 in-memory blocks (glob)
  *Generating BAI (glob)
  *Mapped Reads: 52 (glob)
  *Alignments: 96 (glob)
  *Mapped Bases: 243366 (glob)
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
