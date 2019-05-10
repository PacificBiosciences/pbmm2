pbmm2_cram_test_data = files([
  'baminput.t',
  'biosampleConsensus.t',
  'biosampleSubreads.t',
  'fastxinput.t',
  'idempotence.t',
  'jsoninput.t',
  'misc.t',
  'nohpc.t',
  'overlap.t',
  'readgroup.t',
  'refxml.t',
  'rtc.t',
  'sa.t',
  'splitsample.t',
  'splitsamplejson.t',
  'splitsamplejsonsorted.t',
  'splitsamplesamename.t',
  'splitsamplesorted.t',
  'strip.t',
  'unaligned.t',
  'unrolled.t'
  ])

test(
  'pbmm2_cram test',
  pbmm2_cram_script,
  args : [
    '--xunit-file=cram.xml',
    '--verbose'] +
    pbmm2_cram_test_data,
  env : [
    '__PBTEST_PBMM2_EXE=' + pbmm2_main.full_path(),
    'TESTDIR=' + meson.current_source_dir()],
  timeout : 1800)
