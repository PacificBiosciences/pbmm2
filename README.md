pbmm2
=====

A frontend for Heng Li's speedy minimap2 for PacBio data: native
PacBioBAMs in, native PacBioBAMs out.

Many thanks to Heng Li for a pleasant API experience, and of course to
my colleague Derek Barnett for a similarly pleasant experience (pbbam,
pbcopper).

TODO
----
- Output filter (length and accuracy of alignments)
- Sort and index the output file by default?
- Benchmarks!!

DEPENDENCIES
------------
 - C++14 compiler
 - Meson
 - htslib installed

BUILD
-----

```sh
> mkdir -p build && pushd build && meson && popd
> pushd subprojects/minimap2-2.6 && patch -p1 < ../../m_to_eqx.diff && popd
> pushd build && ninja
> ./src/pbmm2
```
