// Author: Armin Töpfer

#include <functional>
#include <iostream>
#include <memory>
#include <thread>
#include <tuple>
#include <vector>

#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>

#include <pbbam/DataSet.h>
#include <pbbam/DataSetTypes.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>

#include <pbcopper/parallel/WorkQueue.h>

#include <Pbmm2Version.h>

#include "IndexSettings.h"
#include "MM2Helper.h"

#include "IndexWorkflow.h"

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

std::string OutputFilePrefix(const std::string& outputFile)
{
    // Check if output type is a dataset
    const std::string outputExt = Utility::FileExtension(outputFile);
    std::string prefix = outputFile;
    if (outputExt == "xml") {
        boost::ireplace_last(prefix, ".consensusalignmentset.xml", "");
        boost::ireplace_last(prefix, ".alignmentset.xml", "");
    } else if (outputExt == "bam") {
        boost::ireplace_last(prefix, ".bam", "");
        boost::ireplace_last(prefix, ".subreads", "");
    } else {
        PBLOG_FATAL << "Unknown file extension for output file: " << outputFile;
        std::exit(EXIT_FAILURE);
    }
    return prefix;
}
}  // namespace

int IndexWorkflow::Runner(const CLI::Results& options)
{
    std::ofstream logStream;
    {
        const std::string logFile = options["log_file"];
        const Logging::LogLevel logLevel(
            options.IsFromRTC() ? options.LogLevel() : options["log_level"].get<std::string>());

        using Logger = PacBio::Logging::Logger;

        Logger* logger;
        if (!logFile.empty()) {
            logStream.open(logFile);
            logger = &Logger::Default(new Logger(logStream, logLevel));
        } else {
            logger = &Logger::Default(new Logger(std::cout, logLevel));
        }
        PacBio::Logging::InstallSignalHandlers(*logger);
    }

    IndexSettings settings(options);

    if (!CheckPositionalArgs(options.PositionalArguments())) std::exit(EXIT_FAILURE);

    const std::string refFile = options.PositionalArguments()[0];
    const std::string outFile = options.PositionalArguments()[1];

    if (Utility::FileExists(outFile))
        PBLOG_WARN << "Warning: Overwriting existing output file: " << outFile;

    MM2Helper mm2helper(refFile, settings.NumThreads, outFile);

    return EXIT_SUCCESS;
}
}  // namespace minimap2
}  // namespace PacBio
