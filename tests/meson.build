test(
  'pbmm2 formatting check',
  pbmm2_clang_formatter,
  args : [
    '--all'],
  workdir : meson.source_root())

subdir('unit')

pbmm2_cram_test = [
  'baminput',
  'biosampleConsensus',
  'biosampleSubreads',
  'demuxed',
  'empty',
  'enforcedMapping',
  'fastxinput',
  'idempotence',
  'jsoninput',
  'misc',
  'nohpc',
  'overlap',
  'readgroup',
  'refxml',
  'sa',
  'splitsample',
  'splitsamplejson',
  'splitsamplejsonsorted',
  'splitsamplesamename',
  'splitsamplesorted',
  'strip',
  'unaligned',
  'unrolled',
  ]

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

pbmm2_test_samtools = find_program('samtools', required : true)

foreach t : pbmm2_cram_test
  test(
    'pbmm2 cram test - ' + t,
    pbmm2_cram_script,
    args : [
      '--xunit-file=' + join_paths(meson.build_root(), 'pbmm2-cram-' + t.underscorify() + '-test.xml'),
      '--verbose'] + files('cram/' + t + '.t'),
    env : [
      'MINIMAP2=' + pbmm2_test_minimap2_path,
      '__PBTEST_PBMM2_EXE=' + pbmm2_main.full_path(),
      'TESTDIR=' + meson.current_source_dir(),
      'SAMTOOLS=' + pbmm2_test_samtools.path(),
    ],
    timeout : 36000) # with '-O0 -g' tests can be *very* slow
endforeach
