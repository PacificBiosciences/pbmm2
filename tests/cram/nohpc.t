  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" index "$REF" "$CRAMTMP"/index_default.mmi --log-level=DEBUG --preset SUBREAD 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : true (glob)
  $ "$PBMM2" align -j 1 "$IN" "$CRAMTMP"/index_default.mmi "$CRAMTMP"/out_hpc.bam --log-level=DEBUG --preset SUBREAD 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : true (glob)
  $ "$PBMM2" index "$REF" "$CRAMTMP"/index_default.mmi --log-level=DEBUG --preset SUBREAD -u 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : false (glob)
  $ "$PBMM2" align -j 1 "$IN" "$CRAMTMP"/index_default.mmi "$CRAMTMP"/out_hpc.bam --log-level=DEBUG --preset SUBREAD -u 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : false (glob)
  $ "$PBMM2" index "$REF" "$CRAMTMP"/index_default.mmi --log-level=DEBUG --preset ISOSEQ 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : false (glob)
  $ "$PBMM2" align -j 1 "$IN" "$CRAMTMP"/index_default.mmi "$CRAMTMP"/out_hpc.bam --log-level=DEBUG --preset ISOSEQ 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : false (glob)
  $ "$PBMM2" index "$REF" "$CRAMTMP"/index_default.mmi --log-level=DEBUG --preset CCS 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : false (glob)
  $ "$PBMM2" align -j 1 "$IN" "$CRAMTMP"/index_default.mmi "$CRAMTMP"/out_hpc.bam --log-level=DEBUG --preset CCS 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : false (glob)

  $ samtools view -H "$TESTDIR"/data/median.bam > "$CRAMTMP"/single_subread.sam
  $ samtools view -h "$TESTDIR"/data/median.bam | head -n 1 >> "$CRAMTMP"/single_subread.sam
  $ samtools view -bS "$CRAMTMP"/single_subread.sam > "$CRAMTMP"/single_subread.bam
  $ IN="$CRAMTMP"/single_subread.bam

  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/out_hpc.bam --log-level=DEBUG --preset SUBREAD 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : true (glob)
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/out_hpc.bam --log-level=DEBUG --preset SUBREAD -u 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : false (glob)
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/out_hpc.bam --log-level=DEBUG --preset ISOSEQ 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : false (glob)
  $ "$PBMM2" align -j 1 "$IN" "$REF" "$CRAMTMP"/out_hpc.bam --log-level=DEBUG --preset CCS 2>&1 | grep "Homopolymer compressed"
  *Homopolymer compressed : false (glob)
