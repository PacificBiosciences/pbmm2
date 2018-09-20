pbmm2_cram_test_data = files([
  'baminput.t',
  'jsoninput.t',
  'refxml.t'])

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
