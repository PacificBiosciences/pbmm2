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

MM2Helper::MM2Helper(const std::string& refs, const MM2Settings& settings,
                     const std::string& outputMmi)
    : NumThreads{settings.NumThreads}
{
    mm_idxopt_init(&IdxOpts);
    switch (settings.AlignMode) {
        case AlignmentMode::SUBREADS:
            IdxOpts.flag |= MM_I_HPC;
            IdxOpts.k = 19;
            IdxOpts.w = 10;
            break;
        case AlignmentMode::CCS:
            IdxOpts.k = 19;
            IdxOpts.w = 10;
            break;
        case AlignmentMode::ISOSEQ:
            IdxOpts.k = 15;
            IdxOpts.w = 5;
            break;
        default:
            PBLOG_FATAL << "No AlignmentMode --preset selected!";
            std::exit(EXIT_FAILURE);
    }
    if (settings.DisableHPC && IdxOpts.flag & MM_I_HPC) IdxOpts.flag &= ~MM_I_HPC;
    if (settings.Kmer > 0) IdxOpts.k = settings.Kmer;
    if (settings.MinimizerWindowSize > 0) IdxOpts.w = settings.MinimizerWindowSize;
    IdxOpts.batch_size = 0x7fffffffffffffffL;  // always build a uni-part index

    mm_mapopt_init(&MapOpts);
    MapOpts.flag |= MM_F_CIGAR;
    MapOpts.flag |= MM_F_SOFTCLIP;
    MapOpts.flag |= MM_F_LONG_CIGAR;
    MapOpts.flag |= MM_F_EQX;
    MapOpts.flag |= MM_F_NO_PRINT_2ND;
    MapOpts.flag |= MM_F_HARD_MLEVEL;
    MapOpts.mask_level = 0;
    switch (settings.AlignMode) {
        case AlignmentMode::SUBREADS:
            MapOpts.a = 2;
            MapOpts.q = 5;
            MapOpts.q2 = 56;
            MapOpts.e = 4;
            MapOpts.e2 = 1;
            MapOpts.b = 5;
            MapOpts.zdrop = 400;
            MapOpts.zdrop_inv = 50;
            MapOpts.bw = 2000;
            break;
        case AlignmentMode::CCS:
            MapOpts.a = 2;
            MapOpts.q = 5;
            MapOpts.q2 = 56;
            MapOpts.e = 4;
            MapOpts.e2 = 1;
            MapOpts.b = 5;
            MapOpts.zdrop = 400;
            MapOpts.zdrop_inv = 50;
            MapOpts.bw = 2000;
            break;
        case AlignmentMode::ISOSEQ:
            MapOpts.flag |= MM_F_SPLICE | MM_F_SPLICE_FOR | MM_F_SPLICE_FLANK;
            MapOpts.max_gap = 2000;
            MapOpts.max_gap_ref = 200000;
            MapOpts.bw = 200000;
            MapOpts.a = 1;
            MapOpts.b = 2;
            MapOpts.q = 2;
            MapOpts.e = 1;
            MapOpts.q2 = 32;
            MapOpts.e2 = 0;
            MapOpts.zdrop = 200;
            MapOpts.zdrop_inv = 100;
            MapOpts.noncan = 5;
            break;
        default:
            PBLOG_FATAL << "No AlignmentMode --preset selected!";
            std::exit(EXIT_FAILURE);
    }
    if (settings.GapOpenDelete > 0) MapOpts.q = settings.GapOpenDelete;
    if (settings.GapOpenInsert > 0) MapOpts.q2 = settings.GapOpenInsert;
    if (settings.GapExtensionDelete > 0) MapOpts.e = settings.GapExtensionDelete;
    if (settings.GapExtensionInsert > 0) MapOpts.e2 = settings.GapExtensionInsert;
    if (settings.MatchScore > 0) MapOpts.a = settings.MatchScore;
    if (settings.MismatchPenalty > 0) MapOpts.b = settings.MismatchPenalty;
    if (settings.Zdrop > 0) MapOpts.zdrop = settings.Zdrop;
    if (settings.ZdropInv > 0) MapOpts.zdrop_inv = settings.ZdropInv;
    if (settings.NonCanon > 0) MapOpts.noncan = settings.NonCanon;
    if (settings.MaxIntronLength > 0) mm_mapopt_max_intron_len(&MapOpts, settings.MaxIntronLength);
    if (settings.Bandwidth > 0) MapOpts.bw = settings.Bandwidth;
    if (settings.NoSpliceFlank) MapOpts.flag &= ~MM_F_SPLICE_FLANK;

    Idx = std::make_unique<Index>(refs, IdxOpts, NumThreads, outputMmi);
    mm_mapopt_update(&MapOpts, Idx->idx_);
    PBLOG_DEBUG << "Minimap2 parameters";
    PBLOG_DEBUG << "Kmer size              : " << Idx->idx_->k;
    PBLOG_DEBUG << "Minimizer window size  : " << Idx->idx_->w;
    PBLOG_DEBUG << "Homopolymer compressed : " << (Idx->idx_->flag & MM_I_HPC);
    if (outputMmi.empty()) {
        PBLOG_DEBUG << "Deletion gap open      : " << MapOpts.q;
        PBLOG_DEBUG << "Insertion gap open     : " << MapOpts.q2;
        PBLOG_DEBUG << "Deletion gap extension : " << MapOpts.e;
        PBLOG_DEBUG << "Insertion gap extension: " << MapOpts.e2;
        PBLOG_DEBUG << "Match score            : " << MapOpts.a;
        PBLOG_DEBUG << "Mismatch penalty       : " << MapOpts.b;
        PBLOG_DEBUG << "Z-drop                 : " << MapOpts.zdrop;
        PBLOG_DEBUG << "Z-drop inv             : " << MapOpts.zdrop_inv;
        PBLOG_DEBUG << "Bandwidth              : " << MapOpts.bw;
        if (settings.AlignMode == AlignmentMode::ISOSEQ) {
            PBLOG_DEBUG << "Max ref intron length  : " << MapOpts.max_gap_ref;
            PBLOG_DEBUG << "Prefer splice flanks   : " << (!settings.NoSpliceFlank ? "yes" : "no");
        }
    }
}

std::unique_ptr<std::vector<AlignedRecord>> MM2Helper::Align(
    const std::unique_ptr<std::vector<BAM::BamRecord>>& records, const FilterFunc& filter,
    int32_t* alignedReads) const
{
    using namespace PacBio::BAM;
    ThreadBuffer tbuf;
    auto result = std::make_unique<std::vector<AlignedRecord>>();
    result->reserve(records->size());

    for (const auto& record : *records) {
        int numAlns;
        const auto seq = record.Sequence();
        const int qlen = seq.length();
        auto alns = mm_map(Idx->idx_, qlen, seq.c_str(), &numAlns, tbuf.tbuf_, &MapOpts, nullptr);
        bool aligned = false;
        for (int i = 0; i < numAlns; ++i) {
            auto aln = alns[i];
            // if no alignment, continue
            if (aln.p == nullptr) continue;
            // secondary alignment
            if (aln.id != aln.parent) continue;

            aligned = true;
            const int32_t refId = aln.rid;
            const Position refStart = aln.rs;
            const Strand strand = aln.rev ? Strand::REVERSE : Strand::FORWARD;
            const Cigar cigar = RenderCigar(&aln, qlen, MapOpts.flag);
            const uint8_t mapq = aln.mapq;
            auto mapped = BamRecord::Mapped(record, refId, refStart, strand, cigar, mapq);
            mapped.Impl().SetSupplementaryAlignment(aln.sam_pri == 0);
            AlignedRecord alnRec{std::move(mapped)};
            if (filter(alnRec)) result->emplace_back(std::move(alnRec));
        }
        *alignedReads += aligned;
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
    PBLOG_INFO << "Start reading/building index";
    auto rdr =
        mm_idx_reader_open(fname.c_str(), &opts, outputMmi.empty() ? nullptr : outputMmi.c_str());
    if (!rdr) throw std::runtime_error("unable to load reference for indexing!");
    idx_ = mm_idx_reader_read(rdr, numThreads);
    if (!idx_) throw std::runtime_error("unable to index reference!");
    mm_idx_reader_close(rdr);
    PBLOG_INFO << "Finished reading/building index";
}

Index::~Index() { mm_idx_destroy(idx_); }

std::vector<PacBio::BAM::SequenceInfo> Index::SequenceInfos() const
{
    std::vector<PacBio::BAM::SequenceInfo> result;
    for (unsigned i = 0; i < idx_->n_seq; ++i) {
        const std::string name = idx_->seq[i].name;
        const std::string len = std::to_string(idx_->seq[i].len);
        result.emplace_back(
            PacBio::BAM::SequenceInfo(name, len).Checksum("FAKEFAKEFAKEFAKEFAKEFAKEFAKEFAKE"));
    }
    return result;
}

AlignedRecord::AlignedRecord(BAM::BamRecord record) : Record(std::move(record))
{
    ComputeAccuracyBases();
}

void AlignedRecord::ComputeAccuracyBases()
{
    int32_t ins = 0;
    int32_t del = 0;
    int32_t mismatch = 0;
    int32_t n = 0;
    int32_t match = 0;
    for (const auto& cigar : Record.CigarData()) {
        int32_t len = cigar.Length();
        switch (cigar.Type()) {
            case BAM::CigarOperationType::INSERTION:
                ins += len;
                break;
            case BAM::CigarOperationType::DELETION:
                del += len;
                break;
            case BAM::CigarOperationType::SEQUENCE_MISMATCH:
                mismatch += len;
                break;
            case BAM::CigarOperationType::REFERENCE_SKIP:
                break;
            case BAM::CigarOperationType::SEQUENCE_MATCH:
            case BAM::CigarOperationType::ALIGNMENT_MATCH:
                match += len;
                break;
            case BAM::CigarOperationType::PADDING:
            case BAM::CigarOperationType::SOFT_CLIP:
            case BAM::CigarOperationType::HARD_CLIP:
                break;
            case BAM::CigarOperationType::UNKNOWN_OP:
            default:
                PBLOG_FATAL << "UNKNOWN OP";
                std::exit(EXIT_FAILURE);
                break;
        }
    }
    Span = Record.AlignedEnd() - Record.AlignedStart();
    const int32_t nErr = ins + del + mismatch;
    NumAlignedBases = match + ins + mismatch;
    Concordance = 100 * (1.0 - 1.0 * nErr / Span);
}

}  // namespace minimap2
}  // namespace PacBio
