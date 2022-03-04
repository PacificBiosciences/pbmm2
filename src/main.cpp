#include <pbmm2/LibraryInfo.h>
#include "AbortException.h"
#include "AlignSettings.h"
#include "AlignWorkflow.h"
#include "IndexSettings.h"
#include "IndexWorkflow.h"
#include "Pbmm2GlobalVersion.h"

#include <pbcopper/cli2/CLI.h>

#include <iostream>
#include <stdexcept>

PacBio::CLI_v2::MultiToolInterface CreateMultiInterface()
{
    PacBio::CLI_v2::MultiToolInterface mi{"pbmm2", "minimap2 with native PacBio BAM support",
                                          PacBio::Pbmm2::LibraryInfo().Release};

    // clang-format off
    mi.AddTools(
    {
        {"index",
            PacBio::minimap2::IndexSettings::CreateCLI(),
           &PacBio::minimap2::IndexWorkflow::Runner},
        {"align",
            PacBio::minimap2::AlignSettings::CreateCLI(),
           &PacBio::minimap2::AlignWorkflow::Runner}
    });

    mi.HelpFooter(
R"(Typical workflows:
  A. Generate index file for reference and reuse it to align reads
     $ pbmm2 index ref.fasta ref.mmi
     $ pbmm2 align ref.mmi movie.subreads.bam ref.movie.bam

  B. Align reads and sort on-the-fly, with 4 alignment and 2 sort threads
     $ pbmm2 align ref.fasta movie.subreads.bam ref.movie.bam --sort -j 4 -J 2

  C. Align reads, sort on-the-fly, and create PBI
     $ pbmm2 align ref.fasta movie.subreadset.xml ref.movie.alignmentset.xml --sort

  D. Omit output file and stream BAM output to stdout
     $ pbmm2 align hg38.mmi movie1.subreadset.xml | samtools sort > hg38.movie1.sorted.bam

  E. Align CCS fastq input and sort on-the-fly
     $ pbmm2 align ref.fasta movie.Q20.fastq ref.movie.bam --preset CCS --sort --rg '@RG\tID:myid\tSM:mysample')");

    mi.RegisterVersionPrinter(PacBio::Pbmm2::PrintPbmm2VersionMulti);

    // clang-format on
    return mi;
}

int main(int argc, char* argv[]) { return PacBio::CLI_v2::Run(argc, argv, CreateMultiInterface()); }
