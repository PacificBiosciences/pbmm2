  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/sorted.bam --sort
  $ samtools view -H "$CRAMTMP"/sorted.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sorted.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/sorted.bam.bai 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sorted.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/sorted.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/sorted.json 2> /dev/null | wc -l | tr -d ' '
  0
