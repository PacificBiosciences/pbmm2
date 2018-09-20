  $ IN=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/unsorted_header.alignmentset.xml --log-level FATAL
  $ grep -c "PacBio.ReferenceFile.ReferenceFastaFile" $CRAMTMP/unsorted_header.alignmentset.xml
  1
