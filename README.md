<h1 align="center">pbmm2</h1>

A frontend for Heng Li's speedy minimap2 for PacBio data: native
PacBioBAMs in, native PacBioBAMs out.

```sh
> pbmm2 movie.subreads.bam ref.fasta ref.movie.aligned.bam
> pbmm2 movie.subreadset.xml hg38.subreadset.xml hg38.movie.alignmentset.xml
```


## Source Build

### Dependencies
 - C++14 compiler
 - Meson
 - htslib installed

### Build
```sh
> meson build .
> ninja -C build
> ./build/src/pbmm2
```

Collaborators
-------------
Many thanks to Heng Li for a pleasant API experience, of course to
my colleague Derek Barnett for a similarly pleasant experience (pbbam,
pbcopper), and to Lance Hepler for the initial idea and code.