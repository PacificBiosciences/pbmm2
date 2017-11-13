// Copyright (c) 2014-2017, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Armin Töpfer

#include "Settings.h"

namespace PacBio {
namespace minimap2 {
namespace OptionNames {
// clang-format off
const PlainOption NumThreads{
    "numthreads",
    { "j", "num-threads" },
    "Number of Threads",
    "Number of threads to use, 0 means autodetection.",
    CLI::Option::IntType(0)
};
const PlainOption NoPbi{
    "nopbi",
    { "no-pbi" },
    "No PBI",
    "Do not generate a PBI file that is needed for SMRTLink.",
    CLI::Option::BoolType()
};
// clang-format on
}  // namespace OptionNames

Settings::Settings(const PacBio::CLI::Results& options)
    : CLI(options.InputCommandLine())
    , InputFiles(options.PositionalArguments())
    , NoPbi(options[OptionNames::NoPbi])
{
    int requestedNThreads;
    if (options.IsFromRTC()) {
        requestedNThreads = options.NumProcessors();
    } else {
        requestedNThreads = options[OptionNames::NumThreads];
    }
    NumThreads = ThreadCount(requestedNThreads);
}

size_t Settings::ThreadCount(int n)
{
    const int m = std::thread::hardware_concurrency();
    return std::max(1, std::min(m, n));
}

PacBio::CLI::Interface Settings::CreateCLI()
{
    using Option = PacBio::CLI::Option;
    using Task = PacBio::CLI::ToolContract::Task;

    PacBio::CLI::Interface i{"pbmm2", "minimap2 with native PacBio BAM support",
                             "0.0.1"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    i.AddOptions({
        OptionNames::NumThreads,
        OptionNames::NoPbi        
    });

    // clang-format off
    i.AddPositionalArguments({
        { "input", "Source BAM or DATASET", "INPUT" },
        { "reference", "FASTA", "REFERENCE" },
        { "output", "Output BAM", "OUTPUT" }
    });

    const std::string id = "mapping.tasks.pbmm2";
    Task tcTask(id);
    tcTask.NumProcessors(Task::MAX_NPROC);

    tcTask.InputFileTypes({
        {
            "subread_set",
            "SubreadSet",
            "Subread DataSet or .bam file",
            "PacBio.DataSet.SubreadSet"
        },
        {
            "reference_set",
            "ReferenceSet",
            "ReferenceSet or .fasta file",
            "PacBio.DataSet.ReferenceSet"
        }
    });

    tcTask.OutputFileTypes({
        {
            "aligned_bam_output",
            "AlignmentSet",
            "AlignmentSet for output .bam file",
            "PacBio.DataSet.AlignmentSet",
            "pbmm2_output"
        }
    });

    CLI::ToolContract::Config tcConfig(tcTask);
    i.EnableToolContract(tcConfig);
    // clang-format on

    return i;
}
}
}  // ::PacBio::Lima
