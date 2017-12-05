#############
# libraries #
#############

pbmm2_cpp_sources = files([
  'Settings.cpp',
  'Workflow.cpp',])

pbmm2_lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

pbmm2_lib = library(
  'pbmm2',
  pbmm2_cpp_sources,
  install : true,
  soversion : meson.project_version(),
  version : meson.project_version(),
  dependencies : pbmm2_lib_deps,
  include_directories : pbmm2_include_directories,
  cpp_args : pbmm2_warning_flags)


###############
# executables #
###############

# pbmm2
if not meson.is_subproject()
  pbmm2_mains = executable('pbmm2', files(['main.cpp']),
  install : true,
  dependencies : pbmm2_lib_deps,
  include_directories : pbmm2_include_directories,
  link_with : pbmm2_lib,
  cpp_args : pbmm2_warning_flags)
endif