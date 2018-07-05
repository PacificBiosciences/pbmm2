#pragma once

#include <functional>
#include <memory>
#include <vector>

#include <pbbam/BamRecord.h>
#include <pbcopper/logging/Logging.h>

#include <minimap.h>

#include "Index.h"
#include "ThreadBuffer.h"

namespace PacBio {
namespace minimap2 {
namespace {

PacBio::BAM::Cigar RenderCigar(const mm_reg1_t* const r, const int qlen, const int opt_flag)
{
    using PacBio::BAM::Cigar;

    Cigar cigar;

    if (r->p == nullptr) return cigar;

    uint32_t k, clip_len[2];
    clip_len[0] = r->rev ? qlen - r->qe : r->qs;
    clip_len[1] = r->rev ? r->qs : qlen - r->qe;
    const char clip_char = !(opt_flag & MM_F_SOFTCLIP) ? 'H' : 'S'; /* (sam_flag & 0x800) && */

    if (clip_len[0]) cigar.emplace_back(clip_char, clip_len[0]);
    for (k = 0; k < r->p->n_cigar; ++k) {
        cigar.emplace_back("MIDNSHP=XB"[r->p->cigar[k] & 0xf], r -> p -> cigar[k] >> 4);
    }
    if (clip_len[1]) cigar.emplace_back(clip_char, clip_len[1]);

    return cigar;
}
}  // namespace

struct MapOptions
{
    MapOptions()
    {
        mm_mapopt_init(&opts_);
        opts_.flag |= MM_F_CIGAR;
        // opts_.flag |= MM_F_SOFTCLIP;
        opts_.flag |= MM_F_LONG_CIGAR;
        opts_.flag |= MM_F_EQX;
    }

    void Update(const Index& idx) { mm_mapopt_update(&opts_, idx.idx_); }

    mm_mapopt_t opts_;
};

typedef std::unique_ptr<std::vector<PacBio::BAM::BamRecord>> RecordsType;
typedef std::function<bool(const PacBio::BAM::BamRecord&)> FilterFunc;

RecordsType Align(const RecordsType& records, const Index& idx, const MapOptions& mapOpts,
                  const FilterFunc& filter)
{
    using namespace PacBio::BAM;

    ThreadBuffer tbuf;
    auto result = std::make_unique<std::vector<BamRecord>>();
    result->reserve(records->size());

    for (const auto& record : *records) {
        int numAlns;
        const auto seq = record.Sequence();
        const int qlen = seq.length();
        auto alns =
            mm_map(idx.idx_, qlen, seq.c_str(), &numAlns, tbuf.tbuf_, &mapOpts.opts_, nullptr);
        for (int i = 0; i < numAlns; ++i) {
            auto aln = alns[i];
            // if no alignment, continue
            if (aln.p == nullptr) continue;
            const int32_t refId = aln.rid;
            const Position refStart = aln.rs;
            const Strand strand = aln.rev ? Strand::REVERSE : Strand::FORWARD;
            const Cigar cigar = RenderCigar(&aln, qlen, mapOpts.opts_.flag);
            const uint8_t mapq = aln.mapq;
            auto mapped = BamRecord::Mapped(record, refId, refStart, strand, cigar, mapq);
            if (filter(mapped)) {
                result->emplace_back(std::move(mapped));
                break;  // only report the best alignment for now
            }
        }
        // cleanup
        for (int i = 0; i < numAlns; ++i)
            if (alns[i].p) free(alns[i].p);
        free(alns);
    }

    return result;
}
}  // namespace minimap2
}  // namespace PacBio
