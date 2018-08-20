// Author: Armin Töpfer

#include "MM2Helper.h"

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

MM2Helper::MM2Helper(const std::string& refs, const int32_t nthreads, const std::string& outputMmi)
    : NumThreads{nthreads}
{
    mm_idxopt_init(&IdxOpts);
    IdxOpts.flag |= MM_I_HPC;
    IdxOpts.k = 19;
    IdxOpts.w = 10;
    IdxOpts.batch_size = 0x7fffffffffffffffL;  // always build a uni-part index

    mm_mapopt_init(&MapOpts);
    MapOpts.flag |= MM_F_CIGAR;
    MapOpts.flag |= MM_F_SOFTCLIP;
    MapOpts.flag |= MM_F_LONG_CIGAR;
    MapOpts.flag |= MM_F_EQX;
    MapOpts.flag |= MM_F_NO_PRINT_2ND;
    MapOpts.q = 5;
    MapOpts.q2 = 56;
    MapOpts.e = 4;
    MapOpts.e2 = 1;
    MapOpts.b = 5;
    MapOpts.zdrop = 400;
    MapOpts.zdrop_inv = 50;
    MapOpts.bw = 2000;

    Idx = std::make_unique<Index>(refs, IdxOpts, NumThreads, outputMmi);
    mm_mapopt_update(&MapOpts, Idx->idx_);
}

RecordsType MM2Helper::Align(const RecordsType& records, const FilterFunc& filter) const
{
    using namespace PacBio::BAM;

    thread_local ThreadBuffer tbuf;
    auto result = std::make_unique<std::vector<BamRecord>>();
    result->reserve(records->size());

    for (const auto& record : *records) {
        int numAlns;
        const auto seq = record.Sequence();
        const int qlen = seq.length();
        auto alns = mm_map(Idx->idx_, qlen, seq.c_str(), &numAlns, tbuf.tbuf_, &MapOpts, nullptr);
        int numSecondary = 0;
        for (int i = 0; i < numAlns; ++i) {
            auto aln = alns[i];
            // if no alignment, continue
            if (aln.p == nullptr) continue;
            const int32_t refId = aln.rid;
            const Position refStart = aln.rs;
            const Strand strand = aln.rev ? Strand::REVERSE : Strand::FORWARD;
            const Cigar cigar = RenderCigar(&aln, qlen, MapOpts.flag);
            const uint8_t mapq = aln.mapq;
            auto mapped = BamRecord::Mapped(record, refId, refStart, strand, cigar, mapq);
            if (filter(mapped)) {
                if (numSecondary++ < 5)
                    result->emplace_back(std::move(mapped));
                else
                    break;
            }
        }
        // cleanup
        for (int i = 0; i < numAlns; ++i)
            if (alns[i].p) free(alns[i].p);
        free(alns);
    }

    return result;
}

std::vector<PacBio::BAM::SequenceInfo> MM2Helper::SequenceInfos() const
{
    return Idx->SequenceInfos();
}

Index::Index(const std::string& fname, const mm_idxopt_t& opts, const int32_t& numThreads,
             const std::string& outputMmi)
    : idx_{nullptr}
{
    PBLOG_INFO << "Start reading and indexing references";
    auto rdr =
        mm_idx_reader_open(fname.c_str(), &opts, outputMmi.empty() ? nullptr : outputMmi.c_str());
    if (!rdr) throw std::runtime_error("unable to load reference for indexing!");
    idx_ = mm_idx_reader_read(rdr, numThreads);
    if (!idx_) throw std::runtime_error("unable to index reference!");
    mm_idx_reader_close(rdr);
    PBLOG_INFO << "Finished reading and indexing references";
}

Index::~Index() { mm_idx_destroy(idx_); }

std::vector<PacBio::BAM::SequenceInfo> Index::SequenceInfos() const
{
    std::vector<PacBio::BAM::SequenceInfo> result;
    for (unsigned i = 0; i < idx_->n_seq; ++i) {
        const std::string name = idx_->seq[i].name;
        const std::string len = std::to_string(idx_->seq[i].len);
        result.emplace_back(PacBio::BAM::SequenceInfo(name, len));
    }
    return result;
}

}  // namespace minimap2
}  // namespace PacBio
