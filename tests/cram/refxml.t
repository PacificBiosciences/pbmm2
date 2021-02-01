  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align "$IN" "$REF" "$CRAMTMP"/unsorted_header.alignmentset.xml --log-level FATAL
  $ grep -c "PacBio.ReferenceFile.ReferenceFastaFile" "$CRAMTMP"/unsorted_header.alignmentset.xml
  1
