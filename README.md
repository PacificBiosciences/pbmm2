<h1 align="center"><img width="200px" src="img/pbmm2.png"/></h1>
<h1 align="center">pbmm2</h1>
<p align="center">A minimap2 frontend for PacBio data:
native PacBio BAM in ⇨ native PacBio BAM out.</p>

***

## Availability
The latest pre-release, unstable, experts-only linux/mac binaries can be installed via [bioconda](https://bioconda.github.io/).

    conda install pbmm2

These binaries are not ISO compliant.
For research only.
Not for use in diagnostics procedures.

No support for source builds.
No support via mail to developers.
Please *do not* contact a PacBio Field Applications Scientist or PacBio Customer Service for assistance.
Please file GitHub issues for problems and questions!

**This is an early beta!** Expect extreme changes and different output between
versions until release of the first stable release.
Furthermore, the command-line options are not stable yet,
and can change at any point, do not rely on it yet.

## Scope
_pbmm2_ is a wrapper for [minimap2](https://github.com/lh3/minimap2).
Its purpose is to support native PacBio BAM in- and output and provide sets of
recommended parameters. Extensive testing is yet to be performed before _pbmm2_
becomes an officially recommended PacBio aligner; until then, please use BLASR
if you need ISO compliant tools and official PacBio support.

## Usage
_pbmm2_ offers following tools

```
Tools:
    index      Index reference and store as .mmi file
    align      Align subreads or ccs reads to a reference
```

### Index
Indexing is optional, but recommended it you use the same reference multiple times.
```
Usage: pbmm2 index [options] <ref.fa|xml> <out.mmi>
```

### Align
The output argument is optional. If not provided, BAM output is streamed to stdout.
```
Usage: pbmm2 align [options] <in.subreads.bam|xml> <ref.fa|xml|mmi> [out.aligned.bam|xml]
```

## FAQ
### How can I get sorted alignments for polishing?
```sh
> pbmm2 align movie1.subreadset.xml hg38.mmi | samtools sort > hg38.movie1.sorted.bam
> pbindex hg38.movie1.sorted.bam
```
Do not forget to provide `samtools sort` sufficient number of threads and memory
per thread.

### What seeding and alignment parameters does _pbmm2_ use?
The current default mode uses following options:
```
-x map-pb -a --eqx -L -O 5,56 -E 4,1 -B 5 --secondary=no -z 400,50 -r 2k -Y
```

Explained:
* `-x map-pb` provides preset seeding parameters optimized for PacBio reads `-H -k 19 -w 10`
* `-a` produces SAM output, converted to BAM internally
* `--eqx` uses `X`/`=` extended CIGAR strings
* `-L` supports long alignments
* `-O 5,56 -E 4,1 -B 5` approximates the convex gap costs of `ngmlr`
* `--secondary=no` outputs only primary and supplementary alignments
* `-z 400,50` enables alignment of short inversions
* `-r 2k` increases alignment bandwidth to span large insertions and deletions
* `-Y` prevents hard clipping (**required** for `pbsv`)

### Can I override seeding parameters for the index?
If you want to use different seeding parameters `-k` and `-w`, please create
a `.mmi` index with
```sh
> minimap2 -d ref.mmi ref.fasta -H -k X -w Y
```
and provide that
index to `pbmm2 align`.

### Can I override alignment parameters?
No. If you want to experiment with different parameter sets, please use
_minimap2_ in combination with [_pbbamify_](https://github.com/PacificBiosciences/pbbam/wiki/pbbamify) (`conda install pbbam`).
We are open to parameter set suggestions.

### Are more parameter sets planned?
Yes. The current set has been tested for diploid genome mapping.
A prime example for another parameter set is isoform mapping.
This is work in progress.

### Why is the output different from BLASR?
As for any two alignments of the same data with different mappers, alignments
will differ. This is because of many reasons, but mainly a combination of
different scoring functions and seeding techniques.

## Acknowledgements
Many thanks to Heng Li for a pleasant API experience and
to Lance Hepler for the initial idea and code.

## Disclaimer
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
