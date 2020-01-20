pbmm2_cram_test_data = files([
  'baminput.t',
  'biosampleConsensus.t',
  'biosampleSubreads.t',
  'demuxed.t',
  'empty.t',
  'fastxinput.t',
  'idempotence.t',
  'jsoninput.t',
  'misc.t',
  'nohpc.t',
  'overlap.t',
  'readgroup.t',
  'refxml.t',
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

pbmm2_test_minimap2 = find_program('minimap2', required : false)
if pbmm2_test_minimap2.found()
  pbmm2_test_minimap2_path = pbmm2_test_minimap2.path()
else
  message('Using minimap2 from subproject')

  if pbmm2_minimap2_dep.type_name() == 'internal'
    pbmm2_test_minimap2 = pbmm2_minimap2_subproject.get_variable('mm2_main')
    pbmm2_test_minimap2_path = pbmm2_test_minimap2.full_path()
  else
    error('Could not find minimap2 executable, which is required for tests')
  endif
endif

test(
  'pbmm2_cram test',
  pbmm2_cram_script,
  args : [
    '--xunit-file=' + join_paths(meson.build_root(), 'pbmm2-cram-tests.xml'),
    '--verbose'] +
    pbmm2_cram_test_data,
  env : [
    'MINIMAP2=' + pbmm2_test_minimap2_path,
    '__PBTEST_PBMM2_EXE=' + pbmm2_main.full_path(),
    'TESTDIR=' + meson.current_source_dir()],
  timeout : 1800)
