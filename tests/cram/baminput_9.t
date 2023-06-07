  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/sorteddscsi.alignmentset.xml --preset SUBREAD --sort --bam-index CSI 2> "$CRAMTMP"/sorteddscsi.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/sorteddscsi.err
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ samtools view -H "$CRAMTMP"/sorteddscsi.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sorteddscsi.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sorteddscsi.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/sorteddscsi.bam.csi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sorteddscsi.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sorteddscsi.json 2> /dev/null | wc -l | tr -d ' '
  0
  $ grep PacBio.Index.BamIndex "$CRAMTMP"/sorteddscsi.*.xml | wc -l | tr -d ' '
  0
  $ grep PacBio.Index.CsiIndex "$CRAMTMP"/sorteddscsi.*.xml
  *MetaType="PacBio.Index.CsiIndex" ResourceId="*.bam.csi"* (glob)
