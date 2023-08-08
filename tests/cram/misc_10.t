  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset CCS --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: CCS / HiFi (glob)
  * Kmer size              : 19 (glob)
  * Minimizer window size  : 19 (glob)
  * Homopolymer compressed : false (glob)

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset hifi --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: CCS / HiFi (glob)
  * Kmer size              : 19 (glob)
  * Minimizer window size  : 19 (glob)
  * Homopolymer compressed : false (glob)

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset HiFi --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: CCS / HiFi (glob)
  * Kmer size              : 19 (glob)
  * Minimizer window size  : 19 (glob)
  * Homopolymer compressed : false (glob)

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset isoseq --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: ISOSEQ (glob)
  * Kmer size              : 15 (glob)
  * Minimizer window size  : 5 (glob)
  * Homopolymer compressed : false (glob)

  $ "$PBMM2" index "$REF" "$CRAMTMP"/ccs.mmi --preset sUBreAd --log-level DEBUG 2>&1| grep DEBUG
  * Minimap2 parameters based on preset: SUBREAD (glob)
  * Kmer size              : 19 (glob)
  * Minimizer window size  : 10 (glob)
  * Homopolymer compressed : true (glob)
