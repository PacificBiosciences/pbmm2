Basic functionality of option `-N <int>`

  $ "$PBMM2" align -j 1 --preset HiFi \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.141297211.bc2049--bc2049.reference.fasta \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.141297211.bc2049--bc2049.bam \
  > out.bam

  $ samtools view out.bam | cut -f 1-5 | sed 's/\t/ /g'
  m64385e_230506_131841/141297211/ccs 0 chr4 5061 60
  m64385e_230506_131841/141297211/ccs 2048 chr4 8476 60
  m64385e_230506_131841/141297211/ccs 2064 chr16 11579 1
  m64385e_230506_131841/141297211/ccs 2064 chr16 11681 60
  m64385e_230506_131841/141297211/ccs 2064 chr16 12205 14
  m64385e_230506_131841/141297211/ccs 2064 chr16 12060 12
  m64385e_230506_131841/141297211/ccs 2064 chr16 11846 22

  $ "$PBMM2" align -j 1 --preset HiFi -N 1 \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.141297211.bc2049--bc2049.reference.fasta \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.141297211.bc2049--bc2049.bam \
  > out.n1.bam

  $ samtools view out.n1.bam | cut -f 1-5 | sed 's/\t/ /g'
  m64385e_230506_131841/141297211/ccs 0 chr4 5061 60

  $ "$PBMM2" align -j 1 --preset HiFi -N 4 \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.141297211.bc2049--bc2049.reference.fasta \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.141297211.bc2049--bc2049.bam \
  > out.n4.bam

  $ samtools view out.n4.bam | cut -f 1-5 | sed 's/\t/ /g'
  m64385e_230506_131841/141297211/ccs 0 chr4 5061 60
  m64385e_230506_131841/141297211/ccs 2048 chr4 8476 60
  m64385e_230506_131841/141297211/ccs 2064 chr16 11579 1
  m64385e_230506_131841/141297211/ccs 2064 chr16 11681 60
