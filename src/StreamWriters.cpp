// Author: Armin Töpfer

#include <memory>
#include <thread>

#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <pbbam/PbiFile.h>
#include <pbcopper/PbcopperMakeUnique.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/FileUtils.h>
#include <pbcopper/utility/Stopwatch.h>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include "InputOutputUX.h"
#include "Timer.h"
#include "bam_sort.h"

#include "StreamWriters.h"

namespace PacBio {
namespace minimap2 {
namespace {
std::string CreateTmpFile(const std::string& outFile)
{
    std::string pipeName;
    int counter = 0;
    do {
        pipeName = std::tmpnam(nullptr);
        if (++counter == 100) break;
        PBLOG_DEBUG << "Trying pipe " << pipeName;
    } while (Utility::FileExists(pipeName));
    if (counter == 100) {
        pipeName = outFile + ".pbmm2.pipe";
        PBLOG_DEBUG << "Trying pipe " << pipeName;
    }
    return pipeName;
}

void PrintErrorAndAbort(int error)
{
    if (error == EACCES) {
        PBLOG_FATAL << "Pipe error: "
                    << "A component of the path prefix denies search permission, or write "
                       "permission is denied on the parent directory of the FIFO to be "
                       "created.";
    }
    if (error == EEXIST) {
        PBLOG_FATAL << "Pipe error: "
                    << "The named file already exists. Please remove file!";
    }
    if (error == ELOOP) {
        PBLOG_FATAL << "Pipe error: "
                    << "A loop exists in symbolic links encountered during resolution of "
                       "the path argument.";
    }
    if (error == ENAMETOOLONG) {
        PBLOG_FATAL << "Pipe error: "
                    << "The length of the path argument exceeds {PATH_MAX} or a pathname "
                       "component is longer than {NAME_MAX}.";
    }
    if (error == ENOENT) {
        PBLOG_FATAL << "Pipe error: "
                    << "A component of the path prefix specified by path does not name an "
                       "existing directory or path is an empty string.";
    }
    if (error == ENOSPC) {
        PBLOG_FATAL << "Pipe error: "
                    << "The directory that would contain the new file cannot be extended "
                       "or the file system is out of file-allocation resources.";
    }
    if (error == ENOTDIR) {
        PBLOG_FATAL << "Pipe error: "
                    << "A component of the path prefix is not a directory.";
    }
    if (error == EROFS) {
        PBLOG_FATAL << "Pipe error: "
                    << "The named file resides on a read-only file system.";
    }
    if (error == ELOOP) {
        PBLOG_FATAL << "Pipe error: "
                    << "More than {SYMLOOP_MAX} symbolic links were encountered during "
                       "resolution of the path argument.";
    }
    if (error == ENAMETOOLONG) {
        PBLOG_FATAL << "Pipe error: "
                    << "As a result of encountering a symbolic link in resolution of the "
                       "path argument, the length of the substituted pathname string "
                       "exceeded {PATH_MAX}";
    }
    std::exit(EXIT_FAILURE);
}
}  // namespace

StreamWriter::StreamWriter(BAM::BamHeader header, const std::string& outPrefix, bool sort,
                           int sortThreads, int numThreads, int64_t sortMemory,
                           const std::string& sample, const std::string& infix)
    : sort_(sort)
    , sample_(sample)
    , outPrefix_(outPrefix)
    , sortThreads_(sortThreads)
    , numThreads_(numThreads)
    , sortMemory_(sortMemory)
    , header_(std::move(header))
{
    if (outPrefix_ != "-") {
        if (sample_.empty())
            finalOutputPrefix_ = outPrefix_;
        else
            finalOutputPrefix_ = outPrefix_ + '.' + infix;
        finalOutputName_ = finalOutputPrefix_ + ".bam";
    }

    if (!sample_.empty()) {
        auto rgs = header_.ReadGroups();
        header_.ClearReadGroups();
        for (auto& rg : rgs)
            if (rg.Sample() == sample_) header_.AddReadGroup(rg);
    }

    BAM::BamWriter::Config bamWriterConfig;
    std::string outputFile;
    if (sort_) {
        pipeName_ = CreateTmpFile(finalOutputName_);
        int pipe = mkfifo(pipeName_.c_str(), 0666);
        if (pipe == -1) {
            int error = errno;
            PBLOG_FATAL << "Could not open pipe! File name: " << pipeName_;
            PrintErrorAndAbort(error);
        }

        sortThread_ = std::make_unique<std::thread>([&]() {
            int numFiles = 0;
            int numBlocks = 0;
            bam_sort(pipeName_.c_str(), finalOutputName_.c_str(), sortThreads_,
                     sortThreads_ + numThreads_, sortMemory_, &numFiles, &numBlocks);
            PBLOG_INFO << "Merged sorted output from " << numFiles << " files and " << numBlocks
                       << " in-memory blocks";
        });

        outputFile = pipeName_;
        bamWriterConfig.useTempFile = false;
    } else {
        outputFile = finalOutputName_;
        bamWriterConfig.useTempFile = true;
    }
    bamWriter_ = std::make_unique<BAM::BamWriter>(outputFile, header_, bamWriterConfig);
}

void StreamWriter::Write(const BAM::BamRecord& r) const
{
    if (!bamWriter_) {
        PBLOG_FATAL << "Nullpointer BamWriter";
        std::exit(EXIT_FAILURE);
    }
    bamWriter_->Write(r);
}

int64_t StreamWriter::Close()
{
    bamWriter_.reset();
    if (sort_) {
        Utility::Stopwatch sortTime;
        unlink(pipeName_.c_str());
        if (sortThread_) sortThread_->join();
        return sortTime.ElapsedMilliseconds();
    } else {
        return 0;
    }
}

std::string StreamWriter::FinalOutputName() { return finalOutputName_; }
std::string StreamWriter::FinalOutputPrefix() { return finalOutputPrefix_; }

StreamWriters::StreamWriters(BAM::BamHeader& header, const std::string& outPrefix,
                             bool splitBySample, bool sort, int sortThreads, int numThreads,
                             int64_t sortMemory)
    : header_(header.DeepCopy())
    , outPrefix_(outPrefix)
    , splitBySample_(splitBySample)
    , sort_(sort)
    , sortThreads_(sortThreads)
    , numThreads_(numThreads)
    , sortMemory_(sortMemory)
{}

StreamWriter& StreamWriters::at(const std::string& infix, const std::string& sample)
{
    static const std::string unsplit = "unsplit";
    if (!splitBySample_) {
        if (sampleNameToStreamWriter.find(unsplit) == sampleNameToStreamWriter.cend())
            sampleNameToStreamWriter.emplace(
                unsplit, std::make_unique<StreamWriter>(header_.DeepCopy(), outPrefix_, sort_,
                                                        sortThreads_, numThreads_, sortMemory_));

        return *sampleNameToStreamWriter.at(unsplit);
    } else {
        if (sampleNameToStreamWriter.find(sample) == sampleNameToStreamWriter.cend())
            sampleNameToStreamWriter.emplace(
                sample,
                std::make_unique<StreamWriter>(header_.DeepCopy(), outPrefix_, sort_, sortThreads_,
                                               numThreads_, sortMemory_, sample, infix));
        return *sampleNameToStreamWriter.at(sample);
    }
}

std::string StreamWriters::WriteDatasetsJson(const BAM::DataSet& inFile,
                                             const std::string& origOutFile,
                                             const std::string& refFile, const bool isFromXML,
                                             const bool outputIsJson, const Summary& s,
                                             const std::string& outPrefix, const bool splitSample)
{
    std::string pbiTiming;
    Timer pbiTimer;
    std::vector<std::string> xmlNames;
    std::vector<std::string> ids;
    for (auto& sample_sw : sampleNameToStreamWriter) {
        BAM::BamFile validationBam(sample_sw.second->FinalOutputName());
        BAM::PbiFile::CreateFrom(validationBam);

        std::string id;
        const auto xmlName = InputOutputUX::CreateDataSet(inFile, refFile, isFromXML,
                                                          sample_sw.second->FinalOutputPrefix(),
                                                          origOutFile, &id, s.NumAlns, s.Bases);
        xmlNames.emplace_back(xmlName);
        ids.emplace_back(id);
    }
    pbiTiming = pbiTimer.ElapsedTime();

    if (outputIsJson || splitSample) {
        JSON::Json datastore;
        datastore["createdAt"] = BAM::ToIso8601(std::chrono::system_clock::now());
        datastore["updatedAt"] = BAM::ToIso8601(std::chrono::system_clock::now());
        datastore["version"] = "0.2.2";
        std::vector<JSON::Json> files;
        int i = 0;
        for (auto& sample_sw : sampleNameToStreamWriter) {
            JSON::Json datastoreFile;
            datastoreFile["createdAt"] = BAM::ToIso8601(std::chrono::system_clock::now());
            datastoreFile["description"] = "Aligned and sorted reads as BAM";
            std::ifstream file(sample_sw.second->FinalOutputName(),
                               std::ios::binary | std::ios::ate);
            datastoreFile["fileSize"] = static_cast<int>(file.tellg());

            switch (inFile.Type()) {
                case BAM::DataSet::TypeEnum::SUBREAD:
                    datastoreFile["fileTypeId"] = "PacBio.DataSet.AlignmentSet";
                    break;
                case BAM::DataSet::TypeEnum::CONSENSUS_READ:
                    datastoreFile["fileTypeId"] = "PacBio.DataSet.ConsensusAlignmentSet";
                    break;
                case BAM::DataSet::TypeEnum::TRANSCRIPT:
                    datastoreFile["fileTypeId"] = "PacBio.DataSet.TranscriptAlignmentSet";
                    break;
                default:
                    throw std::runtime_error("Unsupported input type");
            }

            datastoreFile["isChunked"] = false;
            datastoreFile["modifiedAt"] = BAM::ToIso8601(std::chrono::system_clock::now());
            datastoreFile["name"] = "Aligned reads";
            datastoreFile["path"] = xmlNames[i];
            datastoreFile["sourceId"] = "mapping.tasks.pbmm2_align-out-1";
            boost::uuids::random_generator gen;
            datastoreFile["uniqueId"] = ids[i];
            ++i;
            files.emplace_back(datastoreFile);
        }
        datastore["files"] = files;
        std::ofstream datastoreStream(outPrefix + ".json");
        datastoreStream << datastore.dump(2);
    }
    return pbiTiming;
}

std::string StreamWriters::Close()
{
    CreateEmptyIfNoOutput();
    int64_t ms = 0;
    for (auto& sample_sw : sampleNameToStreamWriter) {
        ms += sample_sw.second->Close();
    }
    if (!sort_)
        return "";
    else
        return Timer::ElapsedTimeFromSeconds(ms * 1e6);
}

void StreamWriters::CreateEmptyIfNoOutput()
{
    if (sampleNameToStreamWriter.empty()) {
        splitBySample_ = false;
        this->at("", "");
    }
}
}  // namespace minimap2
}  // namespace PacBio
