#include <iostream>
#include <memory>
#include <thread>
#include <tuple>
#include <vector>

#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>

#include <pbbam/DataSet.h>
#include <pbbam/DataSetTypes.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>

#include "Index.h"
#include "Mapping.h"
#include "Settings.h"
#include "WorkQueue.h"
#include "Workflow.h"

using namespace PacBio::BAM;
using namespace PacBio::Parallel;

typedef std::vector<BamRecord> Results;

namespace PacBio {
namespace minimap2 {
void WriteRecords(BamWriter& out, Results&& results)
{
    if (!results.empty())
        out.Write(results[0]);
    /*
    for (const auto& aln : results)
        out.Write(aln);
    */
}

void WriterThread(WorkQueue<Results>& queue, std::unique_ptr<BamWriter> out)
{
    while (queue.ConsumeWith(WriteRecords, std::ref(*out)));
}

std::tuple<std::string, std::string, std::string>
CheckPositionalArgs(const std::vector<std::string>& args)
{
    if (args.size() != 3) {
        std::cerr << "ERROR: Please provide all three arguments: input barcodes output"
                  << std::endl;
        exit(1);
    }

    const auto inputFile = args[0];
    if (!Utility::FileExists(inputFile)) {
        std::cerr << "ERROR: Input data file does not exist: " << inputFile << std::endl;
        exit(1);
    }
    BAM::DataSet dsInput(inputFile);
    switch (dsInput.Type()) {
        case BAM::DataSet::TypeEnum::SUBREAD:
        case BAM::DataSet::TypeEnum::CONSENSUS_READ:
            break;
        case BAM::DataSet::TypeEnum::ALIGNMENT:
        case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
        case BAM::DataSet::TypeEnum::BARCODE:
        case BAM::DataSet::TypeEnum::REFERENCE:
        default:
            std::cerr << "ERROR: Unsupported input data file " << inputFile << "\n";
            std::cerr << "       of type " << BAM::DataSet::TypeToName(dsInput.Type()) << std::endl;
            exit(1);
    }

    const auto& referenceFiles = args[1];
    if (!Utility::FileExists(referenceFiles)) {
        std::cerr << "ERROR: Input reference file does not exist: " << referenceFiles << std::endl;
        exit(1);
    }
    BAM::DataSet dsRef(referenceFiles);
    switch (dsRef.Type()) {
        case BAM::DataSet::TypeEnum::REFERENCE:
            break;
        case BAM::DataSet::TypeEnum::BARCODE:
        case BAM::DataSet::TypeEnum::SUBREAD:
        case BAM::DataSet::TypeEnum::ALIGNMENT:
        case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
        case BAM::DataSet::TypeEnum::CONSENSUS_READ:
        default:
            std::cerr << "ERROR: Unsupported reference input file " << referenceFiles << "\n";
            std::cerr << "       of type " << BAM::DataSet::TypeToName(dsRef.Type()) << std::endl;
            exit(1);
    }
    const auto fastaFiles = dsRef.FastaFiles();
    if (fastaFiles.size() != 1) {
        std::cerr << "Only one reference sequence allowed" << std::endl;
            exit(1);
    }

    return {inputFile, fastaFiles.front(), args[2]};
}

std::unique_ptr<BAM::internal::IQuery> BamQuery(const BAM::DataSet& ds)
{
    const auto filter = BAM::PbiFilter::FromDataSet(ds);
    std::unique_ptr<BAM::internal::IQuery> query(nullptr);
    if (filter.IsEmpty())
        query = std::make_unique<BAM::EntireFileQuery>(ds);
    else
        query = std::make_unique<BAM::PbiFilterQuery>(filter, ds);
    return query;
};

void CreateDataSet(const BAM::DataSet& originalInputDataset, const std::string& outputFile, const Settings& settings) {
    using BAM::DataSet;
    std::string metatype;
    std::string outputType;
    const auto type = originalInputDataset.Type();
    switch (type) {
        case BAM::DataSet::TypeEnum::SUBREAD:
            metatype = "PacBio.AlignmentFile.AlignmentBamFile";
            outputType = "alignmentset";
            break;
        case BAM::DataSet::TypeEnum::CONSENSUS_READ:
            metatype = "PacBio.AlignmentFile.ConsensusAlignmentBamFile";
            outputType = "consensusalignmentset";
            break;
        default:
            throw std::runtime_error("Unsupported input type");
    }
    DataSet ds(type);
    ds.Attribute("xmlns:pbdm") = "http://pacificbiosciences.com/PacBioDataModel.xsd";
    ds.Attribute("xmlns:pbmeta") = "http://pacificbiosciences.com/PacBioCollectionMetadata.xsd";
    ds.Attribute("xmlns:pbpn") = "http://pacificbiosciences.com/PacBioPartNumbers.xsd";
    ds.Attribute("xmlns:pbrk") = "http://pacificbiosciences.com/PacBioReagentKit.xsd";
    ds.Attribute("xmlns:pbsample") = "http://pacificbiosciences.com/PacBioSampleInfo.xsd";
    ds.Attribute("xmlns:pbbase") = "http://pacificbiosciences.com/PacBioBaseDataModel.xsd";

    std::string fileName = outputFile;
    if (fileName.find("/") != std::string::npos) {
        std::vector<std::string> splits;
        boost::split(splits, fileName, boost::is_any_of("/"));
        fileName = splits.back();
    }
    BAM::ExternalResource resource(metatype, fileName + ".bam");

    if (!settings.NoPbi) {
        BAM::FileIndex pbi("PacBio.Index.PacBioIndex", fileName + ".bam.pbi");
        resource.FileIndices().Add(pbi);
    }

    ds.ExternalResources().Add(resource);

    ds.Name(fileName);
    ds.TimeStampedName(fileName + "-" + BAM::CurrentTimestamp());
    std::ofstream dsOut(outputFile + "." + outputType + ".xml");
    ds.SaveToStream(dsOut);
}

std::string OutputFilePrefix(const std::string& outputFile) {
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
        std::cerr << "ERROR: Unknown file extension for output file: " << outputFile << std::endl;
        exit(1);
    }
    return prefix;
}

int Workflow::Runner(const CLI::Results& options)
{
    using std::cref;
    using std::move;
    using std::ref;

    Settings settings(options);

    IndexOptions idxOpts(settings.Kmer, settings.Window, settings.NumThreads);
    MapOptions mapOpts;

    BAM::DataSet qryFile;
    std::string refFile, outFile, alnFile;

    std::tie(qryFile, refFile, outFile) = CheckPositionalArgs(options.PositionalArguments());

    const auto outputFilePrefix = OutputFilePrefix(outFile);
    if (Utility::FileExtension(outFile) == "xml") alnFile = outputFilePrefix + ".bam";
    else alnFile = outFile;

    if (Utility::FileExists(alnFile))
        std::cerr << "Warning: Overwriting existing output file: " << alnFile << std::endl;
    if (alnFile != outFile && Utility::FileExists(outFile))
        std::cerr << "Warning: Overwriting existing output file: " << outFile << std::endl;

    {
        Index idx(refFile, idxOpts);
        mapOpts.Update(idx);

        auto qryRdr = BamQuery(qryFile);

        const auto bamFiles = qryFile.BamFiles();
        auto hdr = bamFiles.front().Header();
        for (size_t i = 1; i < bamFiles.size(); ++i)
            hdr += bamFiles.at(i).Header();

        {
            for (const auto si : idx.SequenceInfos())
                hdr.AddSequence(si);
            hdr.AddProgram(ProgramInfo("pbmm2").Name("pbmm2").Version("0.0.1"));
        }

        // we use at least 1 thread for the BamWriter, subtract it here
        WorkQueue<Results> queue(std::max(1, settings.NumThreads - 1));

        std::unique_ptr<BamWriter> out(new BamWriter(alnFile, hdr));
        std::future<void> writer = std::async(std::launch::async, WriterThread, ref(queue), move(out));

        BamRecord rec;
        while (qryRdr->GetNext(rec)) queue.ProduceWith(&Align, rec, cref(idx), cref(mapOpts));

        queue.Finalize();
        writer.wait();
    }

    /* disabled for now, until we can get sorting in
    if (!settings.NoPbi) {
        BAM::BamFile validationBam(alnFile);
        BAM::PbiFile::CreateFrom(validationBam,
                                    BAM::PbiBuilder::CompressionLevel::CompressionLevel_1,
                                    settings.NumThreads);
    }
    */
    if (Utility::FileExtension(outFile) == "xml")
        CreateDataSet(qryFile, outputFilePrefix, settings);

    return EXIT_SUCCESS;
}
} // namespace minimap2
} // namespace PacBio