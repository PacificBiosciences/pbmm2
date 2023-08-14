Basic functionality of option `-N <int>`

  $ "$PBMM2" align -j 1 --preset HiFi \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.73859681.bc2049--bc2049.reference.fasta \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.73859681.bc2049--bc2049.bam \
  > out.bam

  $ samtools view out.bam | cut -f 1-5 | sed 's/\t/ /g'
  m64385e_230506_131841/73859681/ccs 16 chr12 14240 60
  m64385e_230506_131841/73859681/ccs 2064 chr12 14239 60
  m64385e_230506_131841/73859681/ccs 2048 chr12 14434 60
  m64385e_230506_131841/73859681/ccs 2064 chr12 14946 60
  m64385e_230506_131841/73859681/ccs 2048 chr12 16876 60
  m64385e_230506_131841/73859681/ccs 2048 chr12 14245 60

  $ "$PBMM2" align -j 1 --preset HiFi -N 1 \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.73859681.bc2049--bc2049.reference.fasta \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.73859681.bc2049--bc2049.bam \
  > out.n1.bam

  $ samtools view out.n1.bam | cut -f 1-5 | sed 's/\t/ /g'
  m64385e_230506_131841/73859681/ccs 16 chr12 14240 60

  $ "$PBMM2" align -j 1 --preset HiFi -N 4 \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.73859681.bc2049--bc2049.reference.fasta \
  > "$TESTDIR"/data/r64385e_20230505_214554_B01.73859681.bc2049--bc2049.bam \
  > out.n4.bam

  $ samtools view out.n4.bam | cut -f 1-5 | sed 's/\t/ /g'
  m64385e_230506_131841/73859681/ccs 16 chr12 14240 60
  m64385e_230506_131841/73859681/ccs 2064 chr12 14239 60
  m64385e_230506_131841/73859681/ccs 2048 chr12 14434 60
  m64385e_230506_131841/73859681/ccs 2064 chr12 14946 60
