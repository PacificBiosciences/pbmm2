project(
  'pbmm2',
  ['c', 'cpp'],
  version : '0.12.0',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++11',
    'c_std=gnu11'],
  license : 'BSD-3',
  meson_version : '>= 0.46.0')

# CXXFLAGS
pbmm2_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wno-missing-field-initializers',
  '-Wno-c99-extensions',
  '-Wno-deprecated-declarations' # rely on std::tmpnam()
  ]
    if cpp.has_argument(cflag)
      pbmm2_flags += cflag
    endif
endforeach

pbmm2_c_flags = []
c = meson.get_compiler('c')
foreach cflag: [
  '-Wc++-compat',
  '-Wno-unused-parameter',
  '-Wno-unused-variable',
  '-Wno-sign-compare']
    if c.has_argument(cflag)
      pbmm2_c_flags += cflag
    endif
endforeach

if get_option('sse41')
  if cpp.has_argument('-msse4.1')
    pbmm2_flags += '-msse4.1'
  else
    error('Requested SSE4.1, but the C++ compiler does not support -msse4.1')
  endif

  if c.has_argument('-msse4.1')
    pbmm2_c_flags += '-msse4.1'
  else
    error('Requested SSE4.1, but the C compiler does not support -msse4.1')
  endif
endif

# dependencies

## threads
pbmm2_thread_dep = dependency('threads', required : true)

## boost
pbmm2_boost_dep = dependency('boost', required : true)

## zlib
pbmm2_zlib_dep = dependency('zlib', fallback : ['zlib', 'zlib_dep'])

## pbbam
pbmm2_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])

## pbcopper
pbmm2_pbcopper_dep = dependency('pbcopper', fallback : ['pbcopper', 'pbcopper_dep'])

## minimap2
pbmm2_minimap2_dep = dependency('minimap2', fallback : ['minimap2', 'mm2_dep'])

pbmm2_lib_deps = [
  pbmm2_thread_dep,
  pbmm2_boost_dep,
  pbmm2_zlib_dep,
  pbmm2_pbbam_dep,
  pbmm2_pbcopper_dep,
  pbmm2_minimap2_dep]

# include
subdir('include')

# sources
subdir('src')

# tests
if not meson.is_subproject()
  if get_option('tests')
    pbmm2_cram_script = find_program('cram', required : true)
    subdir('tests')
  endif
endif

# dependency info

if (not meson.is_subproject())
  import('pkgconfig').generate(
    libraries : pbmm2_lib,
    version : meson.project_version(),
    name : 'pbmm2',
    filebase : 'pbmm2',
    description : 'PacBio C++ library wrapping minimap2 for PacBio native data formats')
endif

pbmm2_dep = declare_dependency(
  include_directories : pbmm2_include_directories,
  link_with : pbmm2_lib,
  dependencies : pbmm2_lib_deps,
  version : meson.project_version(),
  compile_args : pbmm2_flags)
