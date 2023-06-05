Repeated matches trimming tags all alignments of a record with `rm:i:1` if there is at least one overlap in query positions

  $ "$PBMM2" align -j 1 --preset HiFi \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.141297211.bc2049--bc2049.reference.fasta \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.141297211.bc2049--bc2049.bam \
  > out.bam

  $ samtools view -c out.bam
  7

  $ samtools view out.bam | grep -c -o "rm:i:1"
  7
