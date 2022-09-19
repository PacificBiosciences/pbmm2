  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta
  $ REF2="ecoliK12_pbi_March2013.fsa"
  $ cp "$REF" "$REF2"

  $ "$PBMM2" index -j 1 "$REF2" "ecoli.mmi"
  $ "$PBMM2" align -j 1 "$IN" "$REF2" "$CRAMTMP"/unsorted_pbi.bam --pbi
  $ samtools view -H "$CRAMTMP"/unsorted_pbi.bam | grep "@HD" | grep "unknown" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/unsorted_pbi.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/unsorted_pbi.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsorted_pbi.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsorted_pbi.*.xml 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/unsorted_pbi.json 2> /dev/null | wc -l | tr -d ' '
  0
