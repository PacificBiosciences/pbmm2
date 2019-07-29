// Author: Armin Töpfer

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <pbbam/BamRecord.h>
#include <pbbam/FastaSequence.h>
#include <pbbam/MD5.h>
#include <pbbam/SequenceInfo.h>
#include <pbcopper/data/MappedRead.h>
#include <pbcopper/data/Read.h>
#include <pbcopper/logging/Logging.h>

#include <pbmm2/MM2Settings.h>

// In file included from ../include/pbmm2/MM2Helper.h:20,
//                  from ../src/MM2Helper.cpp:3:
// minimap.h:78:17: warning: ISO C++ forbids flexible array member 'cigar' [-Wpedantic]
//   uint32_t cigar[];
//                  ^
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#include <minimap.h>
#pragma GCC diagnostic pop

extern "C" void mm_idxopt_init(mm_idxopt_t*);

namespace PacBio {
namespace minimap2 {

class AlignedRecord;
struct AlignedRead;
using FilterFunc = std::function<bool(const AlignedRecord&)>;

struct Index
{
    Index(std::vector<BAM::FastaSequence>&& refs, const mm_idxopt_t& opts);
    // If you use this ctor, you have to be very careful, because refs must live
    // as long as this Index instance.
    Index(const std::vector<BAM::FastaSequence>& refs, const mm_idxopt_t& opts);
    Index(const std::string& fname, const mm_idxopt_t& opts, const int32_t& numThreads,
          const std::string& outputMmi = "");

    ~Index();

    void IndexFrom(const std::vector<BAM::FastaSequence>& refs, const mm_idxopt_t& opts);

    std::vector<PacBio::BAM::SequenceInfo> SequenceInfos() const;

    mm_idx_t* idx_;
    const char** seq_ = nullptr;
    const char** name_ = nullptr;
    const std::vector<BAM::FastaSequence> refs_;
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
    MM2Helper(std::vector<BAM::FastaSequence>&& refs, const MM2Settings& settings);
    // If you use this ctor, you have to be very careful, because refs must live
    // as long as this MM2Helper instance!
    MM2Helper(const std::vector<BAM::FastaSequence>& refs, const MM2Settings& settings);
    MM2Helper(const std::string& refs, const MM2Settings& settings,
              const std::string& outputMmi = "");

public:
    // BamRecord API
    std::unique_ptr<std::vector<AlignedRecord>> Align(
        const std::unique_ptr<std::vector<BAM::BamRecord>>& records, const FilterFunc& filter,
        int32_t* alignedReads) const;

    std::vector<AlignedRecord> Align(const BAM::BamRecord& record) const;
    std::vector<AlignedRecord> Align(const BAM::BamRecord& record, const FilterFunc& filter) const;
    std::vector<AlignedRecord> Align(const BAM::BamRecord& record,
                                     std::unique_ptr<ThreadBuffer>& tbuf) const;
    std::vector<AlignedRecord> Align(const BAM::BamRecord& record, const FilterFunc& filter,
                                     std::unique_ptr<ThreadBuffer>& tbuf) const;

    // Read/MappedRead API
    std::unique_ptr<std::vector<AlignedRead>> Align(
        const std::unique_ptr<std::vector<Data::Read>>& records,
        const std::function<bool(const AlignedRead&)>& filter, int32_t* alignedReads) const;

    std::vector<AlignedRead> Align(const Data::Read& record) const;
    std::vector<AlignedRead> Align(const Data::Read& record,
                                   const std::function<bool(const AlignedRead&)>& filter) const;
    std::vector<AlignedRead> Align(const Data::Read& record,
                                   std::unique_ptr<ThreadBuffer>& tbuf) const;
    std::vector<AlignedRead> Align(const Data::Read& record,
                                   const std::function<bool(const AlignedRead&)>& filter,
                                   std::unique_ptr<ThreadBuffer>& tbuf) const;

    std::vector<PacBio::BAM::SequenceInfo> SequenceInfos() const;

private:
    void PreInit(const MM2Settings& settings, std::string* preset);
    void PostInit(const MM2Settings& settings, const std::string& preset,
                  const bool postAlignParameter);

private:
    // this is the actual weight-lifting alignment function
    template <typename In, typename Out>
    std::vector<Out> AlignImpl(const In& record, const std::function<bool(const Out&)>& filter,
                               std::unique_ptr<ThreadBuffer>& tbuf) const;

private:
    mm_idxopt_t IdxOpts;
    mm_mapopt_t MapOpts;
    const int32_t NumThreads;
    std::unique_ptr<Index> Idx;
    AlignmentMode alnMode_;
    const bool trimRepeatedMatches_;
    const int32_t maxNumAlns_;
};

// BamRecord API
template <typename T>
class AlignedRecordImpl
{
public:
    AlignedRecordImpl(T record);

public:
    T Record;
    int32_t NumAlignedBases = 0;
    int32_t Span = 0;
    double Concordance = 0;
    bool IsAligned;

private:
    void ComputeAccuracyBases();
};

class AlignedRecord : public AlignedRecordImpl<BAM::BamRecord>
{
public:
    AlignedRecord(BAM::BamRecord record);
};

// Read/MappedRead API
// compatibility shim to mostly emulate BamRecord's member functions
struct CompatMappedRead : public Data::MappedRead
{
    CompatMappedRead(Data::MappedRead mr, const int32_t refIdArg)
        : Data::MappedRead{std::move(mr)}, refId{refIdArg}
    {}

    int32_t refId;
    bool supplAln;

    const std::string& Sequence(...) const;
    uint8_t MapQuality() const;
    BAM::RecordType Type() const;
    Data::Position QueryStart() const;
    Data::Position QueryEnd() const;
    int32_t ReferenceId() const;
    bool IsMapped() const;
    CompatMappedRead& Impl();
    const CompatMappedRead& Impl() const;
    const Data::Cigar& CigarData() const;
    void SetSupplementaryAlignment(bool supplAlnArg);
    bool IsSupplementaryAlignment() const;

    bool HasTag(const char*) const { return false; }
    void RemoveTag(const char*) const {}
    template <typename... Args>
    void EditTag(const char*, Args...) const
    {}
    template <typename... Args>
    void AddTag(const char*, Args...) const
    {}
};

struct AlignedRead : public AlignedRecordImpl<CompatMappedRead>
{
public:
    AlignedRead(CompatMappedRead record);

    static CompatMappedRead Mapped(Data::Read record, int32_t refId, Data::Position refStart,
                                   Data::Strand strand, Data::Cigar cigar, uint8_t MapQ);
};

}  // namespace minimap2
}  // namespace PacBio
