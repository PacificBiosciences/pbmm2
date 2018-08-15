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

#include <minimap.h>

extern "C" void mm_idxopt_init(mm_idxopt_t*);

namespace PacBio {
namespace minimap2 {
typedef std::unique_ptr<std::vector<PacBio::BAM::BamRecord>> RecordsType;
typedef std::function<bool(const PacBio::BAM::BamRecord&)> FilterFunc;

struct Index
{
    Index(const std::string& fname, const mm_idxopt_t& opts, const int32_t& numThreads);

    ~Index();

    std::vector<PacBio::BAM::SequenceInfo> SequenceInfos() const;

    mm_idx_t* idx_;
    std::map<std::string, std::string> m5m_;
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
    MM2Helper(const std::string& refs, const int32_t nthreads = 3);

public:
    RecordsType Align(const RecordsType& records, const FilterFunc& filter) const;
    std::vector<PacBio::BAM::SequenceInfo> SequenceInfos() const;

private:
    mm_idxopt_t IdxOpts;
    mm_mapopt_t MapOpts;
    const int32_t NumThreads;
    std::unique_ptr<Index> Idx;
};
}  // namespace minimap2
}  // namespace PacBio