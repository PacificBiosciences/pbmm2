# version info
pbmm2_ver_arr = meson.project_version().split('.')
pbmm2_major_version = pbmm2_ver_arr[0]
pbmm2_minor_version = pbmm2_ver_arr[1]
pbmm2_patch_version = pbmm2_ver_arr[2]

pbmm2_config = configuration_data()
pbmm2_config.set('PBMM2_VERSION', meson.project_version())
pbmm2_config.set('MINIMAP2_VERSION', pbmm2_minimap2_dep.version())

pbmm2_version_sources = [

  # old <-- to be removed
  vcs_tag(
    input : 'Pbmm2GitHash.cpp.in',
    output : 'Pbmm2GitHash.cpp',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@PBMM2_GIT_SHA1@'),
  configure_file(
    input : 'Pbmm2Version.cpp.in',
    output : 'Pbmm2Version.cpp',
    configuration : pbmm2_config),
]

pbmm2_gen_headers = [
  # new
  vcs_tag(
    input : 'LibraryGitHash.h.in',
    output : 'LibraryGitHash.h',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@PBMM2_GIT_SHA1@'),
  configure_file(
    input : 'LibraryVersion.h.in',
    output : 'LibraryVersion.h',
    configuration : pbmm2_config),
]

pbmm2_lib_cpp_sources = files(['MM2Helper.cpp'])
pbmm2_lib_cpp_sources += pbmm2_version_sources
pbmm2_lib_cpp_sources += pbmm2_gen_headers

pbmm2_install = (not meson.is_subproject()) or (get_option('default_library') != 'static')

pbmm2_lib = library(
  'pbmm2',
  pbmm2_lib_cpp_sources,
  install : pbmm2_install,
  soversion : meson.project_version(),
  version : meson.project_version(),
  dependencies : [
    pbmm2_zlib_dep,
    pbmm2_pbbam_dep,
    pbmm2_pbcopper_dep,
    pbmm2_minimap2_dep,
    pbmm2_boost_dep],
  include_directories : pbmm2_include_directories,
  cpp_args : pbmm2_flags)

pbmm2_src_include_directories = include_directories('.', '../third-party')

pbmm2_cpp_sources = files([
  '../third-party/bam_sort.c',
  'AlignSettings.cpp',
  'AlignWorkflow.cpp',
  'IndexSettings.cpp',
  'IndexWorkflow.cpp',
  'InputOutputUX.cpp',
  'LibraryInfo.cpp',
  'SampleNames.cpp',
  'StreamWriters.cpp',
  'Timer.cpp'])

pbmm2_cpp_sources += pbmm2_gen_headers

# executables
## pbmm2
pbmm2_main = executable('pbmm2', files(['main.cpp']) + pbmm2_cpp_sources,
  install : not meson.is_subproject(),
  dependencies : pbmm2_lib_deps,
  include_directories : [pbmm2_include_directories, pbmm2_src_include_directories],
  link_with : pbmm2_lib,
  cpp_args : pbmm2_flags)
