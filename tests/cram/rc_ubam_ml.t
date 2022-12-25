  $ IN_FWD="$TESTDIR"/data/m84011_220902_175841_s1.zmw58919191.fwd.hifi_reads.bam
  $ IN_REV="$TESTDIR"/data/m84011_220902_175841_s1.zmw58919191.rev.hifi_reads.bam
  $ REF="$TESTDIR"/data/REF_58919191.fasta

  $ ${PBMM2} align -j 1 ${IN_FWD} ${REF} fwd.aligned.bam
  $ ${PBMM2} align -j 1 ${IN_REV} ${REF} rev.aligned.bam

  $ samtools view ${IN_FWD} | awk '{ print $1 "\t" $2 "\t" substr($10,1,20) "\t" substr($10,length($10)-19); }' | column -t
  m84011_220902_175841_s1/58919191/ccs  4  AAAAGACTGCGCTTTTGTGG  CAAGGTTTTACCTTTTACAA
  $ samtools view ${IN_REV} | awk '{ print $1 "\t" $2 "\t" substr($10,1,20) "\t" substr($10,length($10)-19); }' | column -t
  m84011_220902_175841_s1/58919191/ccs  20  TTGTAAAAGGTAAAACCTTG  CCACAAAAGCGCAGTCTTTT

  $ samtools view fwd.aligned.bam | awk '{ print $1 "\t" $2 "\t" substr($10,1,20) "\t" substr($10,length($10)-19); }' | column -t
  m84011_220902_175841_s1/58919191/ccs  16  TTGTAAAAGGTAAAACCTTG  CCACAAAAGCGCAGTCTTTT
  $ samtools view rev.aligned.bam | awk '{ print $1 "\t" $2 "\t" substr($10,1,20) "\t" substr($10,length($10)-19); }' | column -t
  m84011_220902_175841_s1/58919191/ccs  16  TTGTAAAAGGTAAAACCTTG  CCACAAAAGCGCAGTCTTTT

  $ samtools view fwd.aligned.bam | tr '\t' '\n' | sort > fwd.aligned.sam
  $ samtools view rev.aligned.bam | tr '\t' '\n' | sort > rev.aligned.sam
  $ diff fwd.aligned.sam rev.aligned.sam
