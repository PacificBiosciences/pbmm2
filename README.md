pbmm2
=====

A frontend for Heng Li's speedy minimap2 for PacBio data: native
PacBioBAMs in, native PacBioBAMs out.

Many thanks to Heng Li for a pleasant API experience, of course to
my colleague Derek Barnett for a similarly pleasant experience (pbbam,
pbcopper), and to Lance Hepler for the initial idea and code.

TODO
----
- Sort and index the output file by default?
- Benchmarks!

DEPENDENCIES
------------
 - C++14 compiler
 - Meson
 - htslib installed

BUILD
-----

```sh
> meson build .
> ninja -C build
> ./build/src/pbmm2
```
