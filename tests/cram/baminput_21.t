  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$TESTDIR"/data/median.transcriptset.xml "$REF" "$CRAMTMP"/out_trans_upper.TranscriptAlignmentSet.XmL || echo $?
  $ "$PBMM2" align -j 1 "$TESTDIR"/data/median.transcriptset.xml "$REF" "$CRAMTMP"/out_json_upper.JsON || echo $?
  $ "$PBMM2" align -j 1 "$TESTDIR"/data/median.transcriptset.xml "$REF" "$CRAMTMP"/out_xml_upper.XML 2> "$CRAMTMP"/out_xml_upper.err  || echo $?
  1
  $ cut -f 8 -d '|' < "$CRAMTMP"/out_xml_upper.err
  *Output is XML, but of unknown type! Please use alignmentset.xml, consensusalignmentset.xml, or transcriptalignmentset.xml* (glob)

  $ "$PBMM2" align "$IN" "$REF" "$CRAMTMP"/sorted.bam --sort -j 2 -J 2 -m 100M --log-level INFO 2>&1| grep INFO
  *Using 2 threads for alignments, 2 threads for sorting, and 200M bytes RAM for sorting. (glob)
  *READ input file: *median.bam* (glob)
  *REF  input file: *ecoliK12_pbi_March2013.fasta* (glob)
  *Start reading/building index (glob)
  *Finished reading/building index (glob)
  *Merged sorted output from 0 files and 1 in-memory blocks (glob)
  *Generating BAI (glob)
  *Mapped Reads: 52 (glob)
  *Alignments: 96 (glob)
  *Mapped Bases: 242437 (glob)
  *Mean Gap-Compressed Sequence Identity* (glob)
  *Max Mapped Read Length* (glob)
  *Mean Mapped Read Length* (glob)
  *Index Build/Read Time: * (glob)
  *Alignment Time: * (glob)
  *Sort Merge Time: * (glob)
  *BAI Generation Time: * (glob)
  *Run Time: * (glob)
  *CPU Time: * (glob)
  *Peak RSS: * (glob)
