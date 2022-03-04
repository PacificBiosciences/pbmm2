#include "IndexWorkflow.h"

#include <pbmm2/MM2Helper.h>
#include <pbmm2/Pbmm2Version.h>
#include "AbortException.h"
#include "IndexSettings.h"

#include <pbbam/DataSet.h>
#include <pbbam/DataSetTypes.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/FileUtils.h>
#include <boost/algorithm/string/predicate.hpp>

#include <string>
#include <vector>

namespace PacBio {
namespace minimap2 {
namespace {
bool CheckPositionalArgs(const std::vector<std::string>& args)
{
    if (args.size() != 2) {
        PBLOG_FATAL << "Please provide both arguments: input output!";
        PBLOG_FATAL << "EXAMPLE: pbmm2 index reference.fasta output.mmi";
        return false;
    }

    const auto inputFile = args[0];
    if (!Utility::FileExists(inputFile)) {
        PBLOG_FATAL << "Input data file does not exist: " << inputFile;
        return false;
    }
    BAM::DataSet dsInput(inputFile);
    switch (dsInput.Type()) {
        case BAM::DataSet::TypeEnum::REFERENCE:
            break;
        case BAM::DataSet::TypeEnum::SUBREAD:
        case BAM::DataSet::TypeEnum::CONSENSUS_READ:
        case BAM::DataSet::TypeEnum::ALIGNMENT:
        case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
        case BAM::DataSet::TypeEnum::BARCODE:
        default:
            PBLOG_FATAL << "Unsupported input data file " << inputFile << " of type "
                        << BAM::DataSet::TypeToName(dsInput.Type());
            return false;
    }

    if (!boost::algorithm::ends_with(args[1], ".mmi")) {
        PBLOG_FATAL << "Output file must end with .mmi: " << inputFile;
        return false;
    }

    return true;
}
}  // namespace

int IndexWorkflow::Runner(const CLI_v2::Results& options)
{
    IndexSettings settings(options);

    if (!CheckPositionalArgs(options.PositionalArguments())) return EXIT_FAILURE;

    BAM::DataSet dsRef(options.PositionalArguments()[0]);
    const auto fastaFiles = dsRef.FastaFiles();
    if (fastaFiles.size() != 1) {
        throw AbortException("Only one reference sequence allowed!");
    }

    const std::string refFile = fastaFiles.front();
    const std::string outFile = options.PositionalArguments()[1];

    if (Utility::FileExists(outFile))
        PBLOG_WARN << "Warning: Overwriting existing output file: " << outFile;

    MM2Helper mm2helper(refFile, settings, outFile);

    return EXIT_SUCCESS;
}
}  // namespace minimap2
}  // namespace PacBio
