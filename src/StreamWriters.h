// Author: Armin Töpfer

#pragma once

#include <map>
#include <memory>
#include <string>
#include <tuple>

#include "BamIndex.h"
#include "InputOutputUX.h"

namespace PacBio {
namespace BAM {
class BamWriter;
class BamHeader;
class DataSet;
}  // namespace BAM
namespace minimap2 {
struct Summary
{
    int32_t NumAlns = 0;
    int64_t Bases = 0;
    double Concordance = 0;
    std::vector<int32_t> Lengths;
};

struct StreamWriter
{
    StreamWriter(BAM::BamHeader header, const std::string& outPrefix, bool sort,
                 const BamIndex bamIdx, int sortThreads, int numThreads, int64_t sortMemory,
                 const std::string& sample = "", const std::string& infix = "");

    void Write(const BAM::BamRecord& r) const;

    std::pair<int64_t, int64_t> Close();

    std::string FinalOutputName();
    std::string FinalOutputPrefix();

private:
    bool sort_;
    BamIndex bamIdx_;
    std::string sample_;
    std::string outPrefix_;
    int sortThreads_;
    int numThreads_;
    int64_t sortMemory_;
    std::string pipeName_;
    std::unique_ptr<BAM::BamWriter> bamWriter_;
    std::unique_ptr<std::thread> sortThread_;
    std::string finalOutputName_{"-"};
    std::string finalOutputPrefix_{"-"};
    BAM::BamHeader header_;
};

struct StreamWriters
{
    StreamWriters(BAM::BamHeader& header, const std::string& outPrefix, bool splitBySample,
                  bool sort, const BamIndex bamIdx, int sortThreads, int numThreads,
                  int64_t sortMemory);

    StreamWriter& at(const std::string& infix, const std::string& sample);

    std::string WriteDatasetsJson(const UserIO& uio, const Summary& s, const bool splitSample);
    std::string ForcePbiOutput();

    std::pair<std::string, std::string> Close();

    void CreateEmptyIfNoOutput();

private:
    BAM::BamHeader header_;
    const std::string& outPrefix_;
    bool splitBySample_;
    bool sort_;
    BamIndex bamIdx_;
    int sortThreads_;
    int numThreads_;
    int64_t sortMemory_;
    std::map<std::string, std::string> movieNameToSampleName;
    std::map<std::string, std::unique_ptr<StreamWriter>> sampleNameToStreamWriter;
};
}  // namespace minimap2
}  // namespace PacBio
