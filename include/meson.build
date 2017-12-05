###########
# headers #
###########

if not meson.is_subproject()
  install_headers(
    files([
      'pacbio/minimap2/Index.h',
      'pacbio/minimap2/Mapping.h',
      'pacbio/minimap2/ThreadBuffer.h']),
    subdir : 'pacbio/minimap2')
endif

pbmm2_include_directories = include_directories('.')
