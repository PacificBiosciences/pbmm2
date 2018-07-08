<h1 align="center">pbmm2</h1>

A minimap2 frontend for PacBio data: 
native PacBio BAM in, native PacBio BAM out.

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

## Acknowledgements
Many thanks to Heng Li for a pleasant API experience and 
to Lance Hepler for the initial idea and code.

## Disclaimer
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.