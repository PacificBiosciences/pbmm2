  $ IN="$TESTDIR"/data/median.bam
  $ REF="$TESTDIR"/data/ecoliK12_pbi_March2013.fasta


  $ samtools view -H $IN > single_rc.sam
  $ samtools view $IN | head -n 1 | cut -f 1 > single_rc.record
  $ echo "20" >> single_rc.record
  $ samtools view $IN | head -n 1 | cut -f 3-9 >> single_rc.record
  $ samtools view $IN | head -n 1 | cut -f 10 | rev | tr "ATGC" "TACG" >> single_rc.record
  $ samtools view $IN | head -n 1 | cut -f 11- >> single_rc.record
  $ cat single_rc.record | tr '\n' '\t'  >> single_rc.sam
  $ samtools view -bS single_rc.sam -o single_rc.bam

  $ samtools view -H $IN > single_native.sam
  $ samtools view $IN | head -n 1 >> single_native.sam
  $ samtools view -bS single_native.sam -o single_native.bam

  $ ${PBMM2} align -j 1 single_rc.bam ${REF} single_rc.aligned.bam
  $ ${PBMM2} align -j 1 single_native.bam ${REF} single_native.aligned.bam

  $ samtools view single_rc.bam | awk '{ print $1 "\t" $2 "\t" substr($10,1,20) "\t" substr($10,length($10)-19); }' | column -t
  m54019_180613_153909/4915325/0_1790  20  TTATATGAAGAAAGTGTGCG  TCTGCGCAAATGCGTTGGTG
  $ samtools view single_native.bam | awk '{ print $1 "\t" $2 "\t" substr($10,1,20) "\t" substr($10,length($10)-19); }' | column -t
  m54019_180613_153909/4915325/0_1790  4  CACCAACGCATTTGCGCAGA  CGCACACTTTCTTCATATAA

  $ samtools view single_rc.aligned.bam | awk '{ print $1 "\t" $2 "\t" substr($10,1,20) "\t" substr($10,length($10)-19); }' | column -t
  m54019_180613_153909/4915325/0_1790  0     CACCAACGCATTTGCGCAGA  CGCACACTTTCTTCATATAA
  m54019_180613_153909/4915325/0_1790  2064  TTATATGAAGAAAGTGTGCG  TCTGCGCAAATGCGTTGGTG
  $ samtools view single_native.aligned.bam | awk '{ print $1 "\t" $2 "\t" substr($10,1,20) "\t" substr($10,length($10)-19); }' | column -t
  m54019_180613_153909/4915325/0_1790  0     CACCAACGCATTTGCGCAGA  CGCACACTTTCTTCATATAA
  m54019_180613_153909/4915325/0_1790  2064  TTATATGAAGAAAGTGTGCG  TCTGCGCAAATGCGTTGGTG
