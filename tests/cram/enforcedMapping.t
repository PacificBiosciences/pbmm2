# A                         D
#   ==========  C             ========== F|SSS|     G   |SSS|   H   |SSS|
#           ====================      ==========    ==========  ==========
#   ==========                ==========
# B xxxxxxxx                  E
#
# Legend:
# A, B, C, D, E, F, G, H - Individual distinct sequences in the final FASTA file.
# x - SNPs
# A and B are identical except for SNPs. SNPs are not introduced in the region which overlaps C.
# D and E are completely identical.
# F, G and H are distinct except for the segmental duplication.
# "|SSS|" - Segmental duplication

  $ READS=$TESTDIR/data/enforced_mapping-reads.fasta
  $ REF=$TESTDIR/data/enforced_mapping-ref.fasta
  $ EXPECTED=$TESTDIR/data/enforced_mapping-expected.csv
  $ MAPPING=$TESTDIR/data/enforced_mapping-read_to_contig.csv

  $ $__PBTEST_PBMM2_EXE align $REF $READS $CRAMTMP/enforced0.bam -j 1 --log-level FATAL --enforced-mapping $MAPPING --preset CCS
  $ samtools view $CRAMTMP/enforced0.bam | cut -f 1,3 | tr '\t' " " | sort > $CRAMTMP/enforced0.csv
  $ diff $CRAMTMP/enforced0.csv $EXPECTED

  $ $__PBTEST_PBMM2_EXE index $REF $CRAMTMP/enforced1.mmi --log-level FATAL --preset CCS
  $ $__PBTEST_PBMM2_EXE align $CRAMTMP/enforced1.mmi $READS $CRAMTMP/enforced1.bam -j 1 --log-level FATAL --enforced-mapping $MAPPING --preset CCS
  $ samtools view $CRAMTMP/enforced1.bam | cut -f 1,3 | tr '\t' " " | sort > $CRAMTMP/enforced1.csv
  $ diff $CRAMTMP/enforced1.csv $EXPECTED
