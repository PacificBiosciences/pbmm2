  $ IN="$TESTDIR"/data/median.ccs.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" align -j 1 "$IN" "$REF" full.bam --preset CCS --sort --short-sa-cigar
  $ "$PBMM2" align -j 1 "$IN" "$REF" strip.bam  --preset CCS --sort --strip --short-sa-cigar

  $ samtools view full.bam | head -n 1 | cut -f 12- | tr '\t' '\n' | sort
  NM:i:1
  RG:Z:d38aac3b
  SA:Z:ecoliK12_pbi_March2013,150081,-,1577S1574M2D,60,0;
  ac* (glob)
  ec* (glob)
  fi* (glob)
  fn* (glob)
  fp* (glob)
  ma* (glob)
  mg* (glob)
  np* (glob)
  ri* (glob)
  rn* (glob)
  rp* (glob)
  rq* (glob)
  sn* (glob)
  zm* (glob)

  $ samtools view strip.bam | head -n 1 | cut -f 12- | tr '\t' '\n' | sort
  NM:i:1
  RG:Z:d38aac3b
  SA:Z:ecoliK12_pbi_March2013,150081,-,1577S1574M2D,60,0;
  ac* (glob)
  ec* (glob)
  ma* (glob)
  mg* (glob)
  np* (glob)
  rq* (glob)
  sn* (glob)
  zm* (glob)
