<h1 align="center"><img width="200px" src="img/pbmm2.png"/></h1>
<h1 align="center">pbmm2</h1>
<p align="center">A minimap2 frontend for PacBio data:
native PacBio BAM in ⇨ native PacBio BAM out.</p>

***

_pbmm2_ is a wrapper for [minimap2](https://github.com/lh3/minimap2).
Its purpose is to support native PacBio BAM in- and output and provide sets of
recommended parameters. Extensive testing is yet to be performed before _pbmm2_
becomes an officially recommended PacBio aligner; until then, please use BLASR
if you need ISO compliant tools and official PacBio support.

**This is an early beta!** Expect extreme changes and different output between
versions until release of the first stable release.
Furthermore, the command-line options are not stable yet,
and can change at any point, do not rely on it yet.

## Availability
Latest version can be installed via bioconda package `pbmm2`.

Please refer to our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda)
for information on Installation, Support, License, Copyright, and Disclaimer.

## Usage
_pbmm2_ offers following tools

```
Tools:
    index      Index reference and store as .mmi file
    align      Align PacBio reads to a reference
```

### Index
Indexing is optional, but recommended it you use the same reference multiple times.
```
Usage: pbmm2 index [options] <ref.fa|xml> <out.mmi>
```

**Notes:**
 - If you use an index file, you can't override parameters `-k` and `-w` in `pbmm2 align`!
 - Minimap2 parameter `-H` (homopolymer-compressed k-mer) is always on.
 - You can also use existing minimap2 `.mmi` files in `pbmm2 align`.

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

### What are parameter sets and how can I override them?
Per default, _pbmm2_ uses recommended parameter sets to simplify the plethora
of possible combinations. For this, we currently offer:

```
  --preset  Set alignment mode:
             - "SUBREAD" -k 19 -w 10 -d 5 -D 4 -i 56 -I 1 -A 2 -B 5 -z 400 -Z 50 -r 2000
             - "ISOSEQ" -k 15 -w 5 -d 2 -i 32 -D 1 -I 0 -A 1 -B 2 -z 200 -Z 100 -C 5 -r 200000 -G 200000
            Default ["SUBREAD"]
```

If you want to override any of the parameters of your chosen set,
please use the respective options:

```
  -k   k-mer size (no larger than 28). [-1]
  -w   Minizer window size. [-1]
  -A   Matching score. [-1]
  -B   Mismatch penalty. [-1]
  -d   Deletion gap open penalty. [-1]
  -i   Insertion gap open penalty. [-1]
  -D   Deletion gap extension penalty. [-1]
  -I   Insertion gap extension penalty. [-1]
  -z   Z-drop score. [-1]
  -Z   Z-drop inversion score. [-1]
  -r   Bandwidth used in chaining and DP-based alignment. [-1]
```

For `ISOSEQ`, you can override additional parameters:

```
  -G                  Max intron length (changes -r). [-1]
  -C                  Cost for a non-canonical GT-AG splicing. [-1]
  --no-splice-flank   Do not prefer splice flanks GT-AG.
  ```

If you have suggestions for our default parameters or ideas for a new
parameter set, please open a GitHub issue!

### What about secondary alignments?
We currently only provide primary and supplementary alignments. If you have an
use-case that absolutely needs secondary alignments, please open a GitHub issue!

### I can't find large SVs!
The `--min-accuracy` option, whereas accuracy is defined as

```
    1.0 - (#Deletions + #Insertions + #Mismatches) / MappedReferenceSpan
```

will remove alignments with more unmapped than mapped bases.
You can deactivate this filter with `--min-accuracy 0`.

### Why is the output different from BLASR?
As for any two alignments of the same data with different mappers, alignments
will differ. This is because of many reasons, but mainly a combination of
different scoring functions and seeding techniques.

## Acknowledgements
Many thanks to Heng Li for a pleasant API experience and
to Lance Hepler for the initial idea and code.

## Disclaimer
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
