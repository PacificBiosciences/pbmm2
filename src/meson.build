# version info
pbmm2_ver_arr = meson.project_version().split('.')
pbmm2_major_version = pbmm2_ver_arr[0]
pbmm2_minor_version = pbmm2_ver_arr[1]
pbmm2_patch_version = pbmm2_ver_arr[2]

pbmm2_config = configuration_data()
pbmm2_config.set('PBMM2_VERSION', meson.project_version())
pbmm2_config.set('PBMM2_VERSION_MAJOR', pbmm2_major_version)
pbmm2_config.set('PBMM2_VERSION_MINOR', pbmm2_minor_version)
pbmm2_config.set('PBMM2_VERSION_PATCH', pbmm2_patch_version)

pbmm2_version_sources = [
  vcs_tag(
    input : 'Pbmm2GitHash.cpp.in',
    output : 'Pbmm2GitHash.cpp',
    command : ['git', 'describe', '--always', '--dirty=*'],
    replace_string : '@PBMM2_GIT_SHA1@'),
  configure_file(
    input : 'Pbmm2Version.cpp.in',
    output : 'Pbmm2Version.cpp',
    configuration : pbmm2_config)]

pbmm2_src_include_directories = include_directories('.')

pbmm2_cpp_sources = files([
  'MM2Helper.cpp',
  'AlignSettings.cpp',
  'AlignWorkflow.cpp',
  'IndexSettings.cpp',
  'IndexWorkflow.cpp'])

pbmm2_cpp_sources += pbmm2_version_sources

pbmm2_lib = static_library(
  'pbmm2',
  pbmm2_cpp_sources,
  install : false,
  dependencies : pbmm2_lib_deps,
  include_directories : pbmm2_src_include_directories,
  cpp_args : pbmm2_flags)


# executables
## pbmm2
if not meson.is_subproject()
  pbmm2_mains = executable('pbmm2', files(['main.cpp']),
  install : true,
  dependencies : pbmm2_lib_deps,
  include_directories : pbmm2_src_include_directories,
  link_with : pbmm2_lib,
  cpp_args : pbmm2_flags)
endif
