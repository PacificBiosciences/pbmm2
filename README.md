pbmm2
=====

A frontend for Heng Li's speedy minimap2 for PacBio data: native
PacBioBAMs in, native PacBioBAMs out.

Many thanks to Heng Li for a pleasant API experience, and of course to
my colleague Derek Barnett for a similarly pleasant experience (pbbam,
pbcopper).

TODO
----
- A proper frontend using pbcopper (Hkw from minimap2)
- Output filter (length and accuracy of alignments)
- Sort and index the output file by default?
- Benchmarks!!

BUILD
-----

```sh
> git submodule update --init --remote --recursive
> pushd third-party/pbbam && patch -p1 < ../pbbam_htslib_pkgconfig.diff && popd
> pushd third-party/minimap2/src && patch -p1 < ../m_to_eqx.diff && popd
> mkdir build && pushd build && cmake -GNinja .. && ninja
```
