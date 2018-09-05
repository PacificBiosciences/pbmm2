// Author: Armin Töpfer

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <pbbam/BamRecord.h>
#include <pbbam/FastaReader.h>
#include <pbbam/MD5.h>
#include <pbbam/SequenceInfo.h>
#include <pbcopper/PbcopperMakeUnique.h>
#include <pbcopper/logging/Logging.h>

#include "MM2Settings.h"

#include <minimap.h>

extern "C" void mm_idxopt_init(mm_idxopt_t*);

namespace PacBio {
namespace minimap2 {

class AlignedRecord;
using FilterFunc = std::function<bool(const AlignedRecord&)>;

struct Index
{
    Index(const std::string& fname, const mm_idxopt_t& opts, const int32_t& numThreads,
          const std::string& outputMmi = "");

    ~Index();

    std::vector<PacBio::BAM::SequenceInfo> SequenceInfos() const;

    mm_idx_t* idx_;
};

struct ThreadBuffer
{
    ThreadBuffer() { tbuf_ = mm_tbuf_init(); }

    ~ThreadBuffer() { mm_tbuf_destroy(tbuf_); }

    mm_tbuf_t* tbuf_;
};

class MM2Helper
{
public:
    MM2Helper(const std::string& refs, const MM2Settings& settings,
              const std::string& outputMmi = "");

public:
    std::unique_ptr<std::vector<AlignedRecord>> Align(
        const std::unique_ptr<std::vector<BAM::BamRecord>>& records, const FilterFunc& filter,
        int32_t* alignedReads) const;
    std::vector<PacBio::BAM::SequenceInfo> SequenceInfos() const;

private:
    mm_idxopt_t IdxOpts;
    mm_mapopt_t MapOpts;
    const int32_t NumThreads;
    std::unique_ptr<Index> Idx;
};

class AlignedRecord
{
public:
    AlignedRecord(BAM::BamRecord record);

public:
    BAM::BamRecord Record;
    int32_t NumAlignedBases;
    int32_t Span;
    double Similarity;

private:
    void ComputeAccuracyBases();
};
}  // namespace minimap2
}  // namespace PacBio
