<h1 align="center">pbmm2</h1>
<p align="center">A minimap2 frontend for PacBio data:
native PacBio BAM in, native PacBio BAM out.</p>

***

## Availability
The latest pre-release, unstable, experts-only linux/mac binaries can be installed via [bioconda](https://bioconda.github.io/).

    conda install pbmm2

These binaries are not ISO compliant.
For research only.
Not for use in diagnostics procedures.

Official support is only provided for official and stable
[SMRT Analysis builds](http://www.pacb.com/products-and-services/analytical-software/)
provided by PacBio.
No support for source builds.
No support via mail to developers.

This is an early beta!
Expect extreme changes and different output between versions until release of
the first stable release.

## Usage examples
Align and create index in memory
```sh
> pbmm2 align movie.subreads.bam hg38.fasta hg38.movie.aligned.bam
> pbmm2 align movie.ccs.bam hg38.fasta hg38.movie.aligned.bam
> pbmm2 align movie.subreadset.xml hg38.referenceset.xml hg38.movie.alignmentset.xml
> pbmm2 align movie.consensusreadset.xml hg38.referenceset.xml hg38.movie.alignmentset.xml
```

Optionally, create and reuse an index
```sh
> pbmm2 index hg38.fasta hg38.mmi
> pbmm2 align movie1.subreadset.xml hg38.mmi hg38.movie1.alignmentset.xml
> pbmm2 align movie2.subreadset.xml hg38.mmi hg38.movie2.alignmentset.xml
```

## Acknowledgements
Many thanks to Heng Li for a pleasant API experience and
to Lance Hepler for the initial idea and code.

## Disclaimer
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
