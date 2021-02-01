  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/sorteddsnoidx.alignmentset.xml --sort --bam-index NONE 2> "$CRAMTMP"/sorteddsnoidx.err || echo $?
  $ cut -f 8 -d '|' < "$CRAMTMP"/sorteddsnoidx.err
  - Input is not a dataset, but output is. Please use dataset input for full SMRT Link compatibility!
  $ "$SAMTOOLS" view -H "$CRAMTMP"/sorteddsnoidx.bam | grep "@HD" | grep "coordinate" | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sorteddsnoidx.bam.pbi 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sorteddsnoidx.bam.bai 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/sorteddsnoidx.bam.csi 2> /dev/null | wc -l | tr -d ' '
  0
  $ ls -alh "$CRAMTMP"/sorteddsnoidx.*.xml 2> /dev/null | wc -l | tr -d ' '
  1
  $ ls -alh "$CRAMTMP"/sorteddsnoidx.json 2> /dev/null | wc -l | tr -d ' '
  0
  $ grep PacBio.Index.BamIndex "$CRAMTMP"/sorteddsnoidx.*.xml | wc -l | tr -d ' '
  0
  $ grep PacBio.Index.CsiIndex "$CRAMTMP"/sorteddsnoidx.*.xml | wc -l | tr -d ' '
  0
