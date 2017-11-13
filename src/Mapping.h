#pragma once

#include <vector>

#include <pbbam/BamRecord.h>

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

    if (r->p == nullptr)
        return cigar;

    uint32_t k, clip_len[2];
    clip_len[0] = r->rev? qlen - r->qe : r->qs;
    clip_len[1] = r->rev? r->qs : qlen - r->qe;
    const char clip_char = !(opt_flag & MM_F_SOFTCLIP) ? 'H' : 'S';  /* (sam_flag & 0x800) && */

    if (clip_len[0]) cigar.emplace_back(clip_char, clip_len[0]);
    for (k = 0; k < r->p->n_cigar; ++k)
        cigar.emplace_back("MIDN=X"[r->p->cigar[k] & 0xf], r->p->cigar[k] >> 4);
    if (clip_len[1]) cigar.emplace_back(clip_char, clip_len[1]);

    return cigar;
}

}

struct MapOptions
{
    MapOptions() {
		mm_mapopt_init(&opts_);
		opts_.flag |= MM_F_CIGAR;  // always perform alignment
        opts_.flag |= MM_F_SOFTCLIP;  // always soft-clip
    }

    void Update(const Index& idx) {
        mm_mapopt_update(&opts_, idx.idx_);
    }

    mm_mapopt_t opts_;
};

std::vector<PacBio::BAM::BamRecord> Align(const PacBio::BAM::BamRecord rec,
                                          const Index& idx, const MapOptions& mapOpts)
{
    using namespace PacBio::BAM;

	static constexpr const uint32_t nCigarMax = 65535;
    static thread_local ThreadBuffer tbuf;

    std::vector<BamRecord> result;

    int nAlns;
    const auto seq = rec.Sequence();
    const int qlen = seq.length();
    auto alns = mm_map(idx.idx_, qlen, seq.c_str(), &nAlns, tbuf.tbuf_, &mapOpts.opts_, nullptr);
    for (int i = 0; i < nAlns; ++i) {
        auto aln = alns[i];
        // if no alignment, continue
        if (aln.p == nullptr) continue;
        if (aln.p->n_cigar > nCigarMax) {
            // TODO(lhepler): log this
            continue;
        }
        const int32_t refId = aln.rid;
        const Position refStart = aln.rs;
        const Strand strand = aln.rev ? Strand::REVERSE : Strand::FORWARD;
        const Cigar cigar = RenderCigar(&aln, qlen, mapOpts.opts_.flag);
        const uint8_t mapq = aln.mapq;
        result.emplace_back(BamRecord::Mapped(rec, refId, refStart, strand, cigar, mapq));
        free(aln.p);
    }
    free(alns);

    return result;
}

}}
