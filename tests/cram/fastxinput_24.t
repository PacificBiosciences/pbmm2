  $ BAM="$TESTDIR"/data/m54075_180905_225130.ccs.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ samtools view "$BAM" | awk '{ print "@"$1"\n"$10"\n+\n"$11 }' > "$CRAMTMP"/m54075_180905_225130.fastq
  $ FASTQ="$CRAMTMP"/m54075_180905_225130.fastq

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
