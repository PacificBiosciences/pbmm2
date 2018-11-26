  $ IN=$TESTDIR/data/median.bam
  $ REF=$TESTDIR/data/ecoliK12_pbi_March2013.fasta

  $ $__PBTEST_PBMM2_EXE align $IN $REF $CRAMTMP/sorted_cli.bam --sort
  $ cd $CRAMTMP
  $ cp $TESTDIR/data/rtc.median.json $REF $IN $TESTDIR/data/median.json .
  $ $__PBTEST_PBMM2_EXE align --resolved-tool-contract rtc.median.json 2>&1 | grep -v Requested
  *Using 6 threads for alignments, 2 threads for sorting, and 8G bytes RAM for sorting.* (glob)
  *READ input file: median.json* (glob)
  *REF  input file: ecoliK12_pbi_March2013.fasta* (glob)
  *Setting to SUBREAD preset* (glob)
  *Start reading/building index* (glob)
  *Finished reading/building index* (glob)
  *Cannot find biosample name for movie name m54019_180613_153909! Will use fallback.* (glob)
  *Merged sorted output from 0 files and 1 in-memory blocks* (glob)
  *Mapped Reads:* (glob)
  *Alignments:* (glob)
  *Mapped Bases:* (glob)
  *Mean Mapped Concordance:* (glob)
  *Max Mapped Read Length:* (glob)
  *Mean Mapped Read Length:* (glob)
  *Index Build/Read Time:* (glob)
  *Alignment Time:* (glob)
  *Sort Merge Time:* (glob)
  *PBI Generation Time:* (glob)
  *Run Time:* (glob)
  *CPU Time:* (glob)
  *Peak RSS:* (glob)

  $ samtools view $CRAMTMP/sorted_cli.bam > $CRAMTMP/sorted_cli.sam
  $ samtools view $CRAMTMP/rtc_out.bam > $CRAMTMP/rtc_out.sam
  $ diff $CRAMTMP/sorted_cli.sam $CRAMTMP/rtc_out.sam

  $ cp $TESTDIR/data/rtc.median.1nproc.json .
  $ $__PBTEST_PBMM2_EXE align --resolved-tool-contract rtc.median.1nproc.json 2>&1 | grep AlignSettings
  *Using 1 threads for alignments, 1 threads for sorting, and 4G bytes RAM for sorting.* (glob)
