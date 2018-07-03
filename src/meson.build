pbmm2_src_include_directories = include_directories('.')

pbmm2_cpp_sources = files([
  'Settings.cpp',
  'Workflow.cpp',])

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