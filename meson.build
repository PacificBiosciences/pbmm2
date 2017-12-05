project(
  'pbmm2',
  ['c', 'cpp'],
  version : '0.2.0',
  default_options : [
    'buildtype=release',
    'default_library=static',
    'warning_level=3',
    'cpp_std=c++14',
    'c_std=gnu11'],
  license : 'BSD-3',
  meson_version : '>= 0.41.1')

############
# CXXFLAGS #
############

pbmm2_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wno-terminate',
  '-Wno-missing-field-initializers']
    if cpp.has_argument(cflag)
      pbmm2_warning_flags += cflag
    endif
endforeach

################
# dependencies #
################

# threads
pbmm2_thread_dep = dependency('threads', required : true)

# boost
pbmm2_boost_dep = dependency('boost', required : true)

# zlib
pbmm2_zlib_dep = dependency('zlib', fallback : ['zlib', 'zlib_dep'], static : true)

# pbbam
pbmm2_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'], static : true)

# pbcopper
pbmm2_pbcopper_dep = dependency('pbcopper', fallback : ['pbcopper', 'pbcopper_dep'], static : true)

# minimap2
pbmm2_minimap2_dep = dependency('minimap2', fallback : ['minimap2', 'mm2_dep'], static : true)

pbmm2_lib_deps = [
  pbmm2_thread_dep,
  pbmm2_boost_dep,
  pbmm2_zlib_dep,
  pbmm2_pbbam_dep,
  pbmm2_pbcopper_dep,
  pbmm2_minimap2_dep]

###########
# headers #
###########

subdir('include')

#####################
# sources + library #
#####################

subdir('src')

###################
# dependency info #
###################

if not meson.is_subproject()
  import('pkgconfig').generate(
    libraries : pbmm2_lib,
    version : meson.project_version(),
    name : 'pbmm2',
    filebase : 'pbmm2',
    description : 'minimap2 with native PacBio BAM support')
endif

pbmm2_dep = declare_dependency(
  include_directories : pbmm2_include_directories,
  link_with : pbmm2_lib,
  dependencies : pbmm2_lib_deps,
  version : meson.project_version(),
  compile_args : pbmm2_warning_flags)
