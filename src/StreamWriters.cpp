#include "StreamWriters.h"

#include "AbortException.h"
#include "SampleNames.h"
#include "bam_sort.h"

#include <htslib/sam.h>
#include <pbbam/PbiFile.h>
#include <pbcopper/logging/Logging.h>
#include <pbcopper/utility/FileUtils.h>
#include <pbcopper/utility/Stopwatch.h>
#include <boost/algorithm/string.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>

#include <limits.h>
#include <unistd.h>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <thread>

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
    const char* errMsg = nullptr;
    if (error == EACCES) {
        errMsg =
            "Pipe error: A component of the path prefix denies search permission, or write "
            "permission is denied on the parent directory of the FIFO to be "
            "created.";
    }
    if (error == EEXIST) {
        errMsg =
            "Pipe error: "
            "The named file already exists. Please remove file!";
    }
    if (error == ELOOP) {
        errMsg =
            "Pipe error: "
            "A loop exists in symbolic links encountered during resolution of "
            "the path argument.";
    }
    if (error == ENAMETOOLONG) {
        errMsg =
            "Pipe error: "
            "The length of the path argument exceeds {PATH_MAX} or a pathname "
            "component is longer than {NAME_MAX}.";
    }
    if (error == ENOENT) {
        errMsg =
            "Pipe error: "
            "A component of the path prefix specified by path does not name an "
            "existing directory or path is an empty string.";
    }
    if (error == ENOSPC) {
        errMsg =
            "Pipe error: "
            "The directory that would contain the new file cannot be extended "
            "or the file system is out of file-allocation resources.";
    }
    if (error == ENOTDIR) {
        errMsg =
            "Pipe error: "
            "A component of the path prefix is not a directory.";
    }
    if (error == EROFS) {
        errMsg =
            "Pipe error: "
            "The named file resides on a read-only file system.";
    }
    if (error == ELOOP) {
        errMsg =
            "Pipe error: "
            "More than {SYMLOOP_MAX} symbolic links were encountered during "
            "resolution of the path argument.";
    }
    if (error == ENAMETOOLONG) {
        errMsg =
            "Pipe error: "
            "As a result of encountering a symbolic link in resolution of the "
            "path argument, the length of the substituted pathname string "
            "exceeded {PATH_MAX}";
    }
    throw AbortException(errMsg);
}
}  // namespace

StreamWriter::StreamWriter(BAM::BamHeader header, const std::string& outPrefix, bool sort,
                           const BamIndex bamIdx, int sortThreads, int numThreads,
                           int64_t sortMemory, const std::string& sample, const std::string& infix)
    : sort_(sort)
    , bamIdx_(bamIdx)
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
        for (auto& rg : rgs) {
            if (SampleNames::SanitizeSampleName(rg.Sample()) == sample_) {
                header_.AddReadGroup(rg);
            }
        }
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
            std::string tmpdir;
            if (const char* env_p = std::getenv("TMPDIR")) {
                tmpdir = env_p;
                if (!boost::ends_with(tmpdir, "/")) tmpdir += '/';
            }
            const bool useTmpDir = Utility::FileExists(tmpdir);
            if (useTmpDir) {
                PBLOG_DEBUG << "[TMPDIR] Using directory for sorting: " << tmpdir;
            } else if (tmpdir.empty()) {
                char cwd[PATH_MAX];
                if (getcwd(cwd, sizeof(cwd)) != NULL) {
                    PBLOG_DEBUG
                        << "[TMPDIR] Could not find environment variable, using working directory: "
                        << cwd;
                } else {
                    PBLOG_DEBUG
                        << "[TMPDIR] Could not find environment variable, using working directory";
                }
            } else {
                PBLOG_DEBUG << "[TMPDIR] Specified directory does not exist: " << tmpdir;
            }
            int ret = bam_sort(pipeName_.c_str(), finalOutputName_.c_str(), tmpdir.c_str(),
                               useTmpDir, sortThreads_, sortThreads_ + numThreads_, sortMemory_,
                               &numFiles, &numBlocks);
            if (ret == EXIT_FAILURE) {
                throw AbortException("Fatal error in bam sort. Aborting.");
            }
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
        throw AbortException("Nullpointer BamWriter");
    }
    bamWriter_->Write(r);
}

std::pair<int64_t, int64_t> StreamWriter::Close()
{
    bamWriter_.reset();
    if (sort_) {
        int idxMs = 0;
        Utility::Stopwatch sortTime;
        unlink(pipeName_.c_str());
        if (sortThread_) {
            sortThread_->join();
            if (bamIdx_ != BamIndex::_from_string("NONE") && finalOutputName_ != "-") {
                Utility::Stopwatch time;
                int shift = 0;
                std::string idxType = bamIdx_._to_string();
                switch (bamIdx_) {
                    case BamIndex::BAI:
                        shift = 0;
                        break;
                    case BamIndex::CSI:
                        shift = 14;
                        break;
                    default:
                        throw AbortException("Unexpected index type. Falling back to generate BAI");
                        break;
                }
                PBLOG_INFO << "Generating " << idxType;
                int ret = sam_index_build3(finalOutputName_.c_str(), NULL, shift,
                                           std::min(sortThreads_ + numThreads_, 8));

                switch (ret) {
                    case 0:
                        break;
                    case -2: {
                        std::ostringstream os;
                        os << idxType << " Index Generation: Failed to open file "
                           << finalOutputName_;
                        throw AbortException(os.str());
                    }
                    case -3: {
                        std::ostringstream os;
                        os << idxType
                           << " Index Generation: File is in a format that cannot be "
                              "usefully indexed "
                           << finalOutputName_;
                        throw AbortException(os.str());
                    }
                    case -4: {
                        std::ostringstream os;
                        os << idxType << " Index Generation: Failed to create or write index "
                           << finalOutputName_ + '.' << boost::to_lower_copy(idxType);
                        throw AbortException(os.str());
                    }
                    default: {
                        std::ostringstream os;
                        os << idxType << " Index Generation: Failed to create index for "
                           << finalOutputName_;
                        throw AbortException(os.str());
                    }
                }
                idxMs = time.ElapsedNanoseconds();
            }
        }
        return {sortTime.ElapsedNanoseconds(), idxMs};
    } else {
        return {0, 0};
    }
}

std::string StreamWriter::FinalOutputName() { return finalOutputName_; }
std::string StreamWriter::FinalOutputPrefix() { return finalOutputPrefix_; }

StreamWriters::StreamWriters(BAM::BamHeader& header, const std::string& outPrefix,
                             bool splitBySample, bool sort, const BamIndex bamIdx, int sortThreads,
                             int numThreads, int64_t sortMemory)
    : header_(header.DeepCopy())
    , outPrefix_(outPrefix)
    , splitBySample_(splitBySample)
    , sort_(sort)
    , bamIdx_(bamIdx)
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
                unsplit,
                std::make_unique<StreamWriter>(header_.DeepCopy(), outPrefix_, sort_, bamIdx_,
                                               sortThreads_, numThreads_, sortMemory_));

        return *sampleNameToStreamWriter.at(unsplit);
    } else {
        if (sampleNameToStreamWriter.find(sample) == sampleNameToStreamWriter.cend())
            sampleNameToStreamWriter.emplace(
                sample, std::make_unique<StreamWriter>(header_.DeepCopy(), outPrefix_, sort_,
                                                       bamIdx_, sortThreads_, numThreads_,
                                                       sortMemory_, sample, infix));
        return *sampleNameToStreamWriter.at(sample);
    }
}

std::string StreamWriters::WriteDatasetsJson(const UserIO& uio, const Summary& s,
                                             const bool splitSample)
{
    BAM::DataSet ds;
    try {
        if (uio.isFromJson)
            ds = BAM::DataSet{uio.unpackedFromJson};
        else if (!uio.isFastaInput && !uio.isFastqInput)
            ds = BAM::DataSet{uio.inFile};
    } catch (std::runtime_error& e) {
        throw AbortException(e.what());
    }
    std::string pbiTiming;
    Utility::Stopwatch pbiTimer;
    std::vector<std::string> xmlNames;
    std::vector<std::string> ids;
    for (auto& sample_sw : sampleNameToStreamWriter) {
        BAM::BamFile validationBam(sample_sw.second->FinalOutputName());
        BAM::PbiFile::CreateFrom(validationBam);

        std::string id;
        const auto xmlName = InputOutputUX::CreateDataSet(
            ds, uio.refFile, uio.isFromXML, sample_sw.second->FinalOutputPrefix(), uio.outFile, &id,
            s.NumAlns, s.Bases, uio.bamIndex);
        xmlNames.emplace_back(xmlName);
        ids.emplace_back(id);
    }
    pbiTiming = pbiTimer.ElapsedTime();

    if (uio.isToJson || splitSample) {
        JSON::Json datastore;
        const auto now = BAM::ToIso8601(std::chrono::system_clock::now());
        datastore["createdAt"] = now;
        datastore["updatedAt"] = now;
        datastore["version"] = "0.2.2";
        std::vector<JSON::Json> files;
        int i = 0;
        for (auto& sample_sw : sampleNameToStreamWriter) {
            JSON::Json datastoreFile;
            datastoreFile["createdAt"] = now;
            datastoreFile["description"] = "Aligned and sorted reads as BAM";
            std::ifstream file(sample_sw.second->FinalOutputName(),
                               std::ios::binary | std::ios::ate);
            datastoreFile["fileSize"] = static_cast<int>(file.tellg());

            switch (ds.Type()) {
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
            datastoreFile["modifiedAt"] = now;
            datastoreFile["name"] = "Aligned reads";
            datastoreFile["path"] = xmlNames[i];
            datastoreFile["sourceId"] = "mapping.tasks.pbmm2_align-out-1";
            datastoreFile["uniqueId"] = ids[i];
            ++i;
            files.emplace_back(datastoreFile);
        }
        datastore["files"] = files;
        std::ofstream datastoreStream(uio.outPrefix + ".json");
        datastoreStream << datastore.dump(2);
    }
    return pbiTiming;
}

std::string StreamWriters::ForcePbiOutput()
{
    Utility::Stopwatch pbiTimer;
    for (auto& sample_sw : sampleNameToStreamWriter) {
        BAM::BamFile validationBam(sample_sw.second->FinalOutputName());
        BAM::PbiFile::CreateFrom(validationBam);
    }
    return pbiTimer.ElapsedTime();
}

std::pair<std::string, std::string> StreamWriters::Close()
{
    CreateEmptyIfNoOutput();
    int64_t sortMs = 0;
    int64_t idxMs = 0;
    for (auto& sample_sw : sampleNameToStreamWriter) {
        const auto sort_bai = sample_sw.second->Close();
        sortMs += sort_bai.first - sort_bai.second;
        idxMs += sort_bai.second;
    }
    constexpr int64_t i641 = 1;
    sortMs = std::max(sortMs, i641);
    idxMs = std::max(idxMs, i641);
    if (!sort_)
        return {"", ""};
    else
        return {Utility::Stopwatch::PrettyPrintNanoseconds(sortMs),
                Utility::Stopwatch::PrettyPrintNanoseconds(idxMs)};
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
