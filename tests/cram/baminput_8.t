  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/sortedds.alignmentset.xml --sort 2> "$CRAMTMP"/sortedds.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/sortedds.err
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ "$SAMTOOLS" view -H "$CRAMTMP"/sortedds.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sortedds.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sortedds.bam.bai 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sortedds.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/sortedds.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sortedds.json 2> /dev/null | wc -l | tr -d ' '
  0
  $ grep PacBio.Index.PacBioIndex "$CRAMTMP"/sortedds.*.xml
  *MetaType="PacBio.Index.PacBioIndex" ResourceId="*.bam.pbi"* (glob)
  $ grep PacBio.Index.BamIndex "$CRAMTMP"/sortedds.*.xml
  *MetaType="PacBio.Index.BamIndex" ResourceId="*.bam.bai"* (glob)
  $ grep PacBio.Index.CsiIndex "$CRAMTMP"/sortedds.*.xml | wc -l | tr -d ' '
  0
