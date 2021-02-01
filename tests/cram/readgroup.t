  $ BAM="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$SAMTOOLS" view "$BAM" | awk '{ print ">"$1"\n"$10 }' > "$CRAMTMP"/median.fasta
  $ FASTA="$CRAMTMP"/median.fasta
  $ "$PBMM2" index "$REF" "$CRAMTMP"/ecoli.mmi --log-level FATAL
  $ REF="$CRAMTMP"/ecoli.mmi

  $ "$PBMM2" align "$REF" "$FASTA" "$CRAMTMP"/rg.bam --log-level FATAL
  $ "$SAMTOOLS" view -H "$CRAMTMP"/rg.bam | grep "^@RG"
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)
  $ rm "$CRAMTMP"/rg.bam

  $ "$PBMM2" align "$REF" "$FASTA" "$CRAMTMP"/rg.bam --log-level FATAL --rg '@RG'
  *Invalid @RG line. Missing ID field. Please provide following format: '@RG\tID:xyz\tSM:abc'* (glob)
  [1]

  $ "$PBMM2" align "$REF" "$FASTA" "$CRAMTMP"/rg.bam --log-level FATAL --rg '@XY'
  *Invalid @RG line. Missing ID field. Please provide following format: '@RG\tID:xyz\tSM:abc'* (glob)
  [1]

  $ "$PBMM2" align "$REF" "$FASTA" "$CRAMTMP"/rg.bam --log-level FATAL --sample newSampleName
  $ "$SAMTOOLS" view -H "$CRAMTMP"/rg.bam | grep "^@RG"
  @RG\tID:default\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:newSampleName\tPM:SEQUEL (esc)
  $ rm "$CRAMTMP"/rg.bam

  $ "$PBMM2" align "$REF" "$FASTA" "$CRAMTMP"/rg.bam --log-level FATAL --rg '@RG\tID:myid'
  $ "$SAMTOOLS" view -H "$CRAMTMP"/rg.bam | grep "^@RG"
  @RG\tID:myid\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:UnnamedSample\tPM:SEQUEL (esc)
  $ rm "$CRAMTMP"/rg.bam

  $ "$PBMM2" align "$REF" "$FASTA" "$CRAMTMP"/rg.bam --log-level FATAL --rg '@RG\tID:myid\tSM:mysample'
  $ "$SAMTOOLS" view -H "$CRAMTMP"/rg.bam | grep "^@RG"
  @RG\tID:myid\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:mysample\tPM:SEQUEL (esc)
  $ rm "$CRAMTMP"/rg.bam

  $ "$PBMM2" align "$REF" "$FASTA" "$CRAMTMP"/rg.bam --log-level FATAL --rg '@RG\tID:myid\tSM:mysample' --sample newSampleName
  $ "$SAMTOOLS" view -H "$CRAMTMP"/rg.bam | grep "^@RG"
  @RG\tID:myid\tPL:PACBIO\tDS:READTYPE=UNKNOWN\tPU:default\tSM:newSampleName\tPM:SEQUEL (esc)
  $ rm "$CRAMTMP"/rg.bam

  $ "$PBMM2" align "$REF" $BAM "$CRAMTMP"/rg.bam --log-level FATAL --rg '@RG\tID:myid'
  *Cannot override read groups with BAM input. Remove option --rg.* (glob)
  [1]
