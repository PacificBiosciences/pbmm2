// Author: Armin Töpfer

#include <pbmm2/MM2Helper.h>

#include <iostream>
#include <memory>

#include <pbcopper/data/Cigar.h>
#include <pbcopper/data/Position.h>
#include <pbcopper/data/Strand.h>

namespace PacBio {
namespace minimap2 {
namespace {

Data::Cigar RenderCigar(const mm_reg1_t* const r, const int qlen, const int opt_flag)
{
    using Data::Cigar;

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

Data::Cigar RenderCigar(const mm_reg1_t* const r, const int qlen, const int opt_flag, int newQs,
                        int newQe, int* refStartOffset)
{
    using Data::Cigar;

    Cigar cigar;

    if (r->p == nullptr) return cigar;

    uint32_t k, clip_len[2];
    int origQs = r->qs;
    int origQe = r->qe;
    if (r->rev) {
        std::swap(newQs, newQe);
        newQs = qlen - newQs;
        newQe = qlen - newQe;
        std::swap(origQs, origQe);
        origQs = qlen - origQs;
        origQe = qlen - origQe;
    }
    clip_len[0] = newQs;
    clip_len[1] = qlen - newQe;
    const char clip_char = !(opt_flag & MM_F_SOFTCLIP) ? 'H' : 'S'; /* (sam_flag & 0x800) && */

    if (clip_len[0]) cigar.emplace_back(clip_char, clip_len[0]);
    int position = origQs;
    int refSpace = 0;
    for (k = 0; k < r->p->n_cigar; ++k) {
        const char cigarChar = "MIDNSHP=XB"[r->p->cigar[k] & 0xf];
        int used = 0;
        for (size_t l = 0; l<r->p->cigar[k]>> 4; ++l) {
            switch (cigarChar) {
                case 'M':
                case '=':
                case 'X':
                    if ((position < newQs)) ++refSpace;
                    /* Falls through. */
                case 'I':
                    if (position >= newQs && position < newQe) ++used;
                    ++position;
                    break;
                case 'D':
                case 'N':
                    if ((position < newQs)) ++refSpace;
                    if (position >= newQs && position < newQe) ++used;
                    break;
                case 'S':
                case 'H':
                    PBLOG_FATAL << "Cigar should not occur " << cigarChar;
                    std::exit(EXIT_FAILURE);
                default:
                    PBLOG_FATAL << "Unknown cigar " << cigarChar;
                    std::exit(EXIT_FAILURE);
                    break;
            }
        }
        if (used > 0) cigar.emplace_back(cigarChar, used);
    }
    *refStartOffset = refSpace;
    if (clip_len[1]) cigar.emplace_back(clip_char, clip_len[1]);

    return cigar;
}
}  // namespace

MM2Helper::MM2Helper(const std::string& refs, const MM2Settings& settings,
                     const std::string& outputMmi)
    : NumThreads{settings.NumThreads}
    , alnMode_(settings.AlignMode)
    , trimRepeatedMatches_(!settings.NoTrimming)
    , maxNumAlns_(settings.MaxNumAlns)
{
    std::string preset;
    PreInit(settings, &preset);
    Idx = std::make_unique<Index>(refs, IdxOpts, NumThreads, outputMmi);
    PostInit(settings, preset, outputMmi.empty());
}
MM2Helper::MM2Helper(const std::vector<BAM::FastaSequence>& refs, const MM2Settings& settings)
    : NumThreads{settings.NumThreads}
    , alnMode_(settings.AlignMode)
    , trimRepeatedMatches_(!settings.NoTrimming)
    , maxNumAlns_(settings.MaxNumAlns)
{
    std::string preset;
    PreInit(settings, &preset);
    Idx = std::make_unique<Index>(refs, IdxOpts);
    PostInit(settings, preset, true);
}
MM2Helper::MM2Helper(std::vector<BAM::FastaSequence>&& refs, const MM2Settings& settings)
    : NumThreads{settings.NumThreads}
    , alnMode_(settings.AlignMode)
    , trimRepeatedMatches_(!settings.NoTrimming)
    , maxNumAlns_(settings.MaxNumAlns)
{
    std::string preset;
    PreInit(settings, &preset);
    Idx = std::make_unique<Index>(std::move(refs), IdxOpts);
    PostInit(settings, preset, true);
}
void MM2Helper::PreInit(const MM2Settings& settings, std::string* preset)
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
        case AlignmentMode::UNROLLED:
            IdxOpts.flag |= MM_I_HPC;
            IdxOpts.k = 15;
            IdxOpts.w = 15;
            break;
        default:
            PBLOG_FATAL << "No AlignmentMode --preset selected!";
            std::exit(EXIT_FAILURE);
    }
    if (settings.DisableHPC && IdxOpts.flag & MM_I_HPC) IdxOpts.flag &= ~MM_I_HPC;
    if (settings.Kmer >= 0) IdxOpts.k = settings.Kmer;
    if (settings.MinimizerWindowSize >= 0) IdxOpts.w = settings.MinimizerWindowSize;
    IdxOpts.batch_size = 0x7fffffffffffffffL;  // always build a uni-part index

    mm_mapopt_init(&MapOpts);
    MapOpts.flag |= MM_F_CIGAR;
    MapOpts.flag |= MM_F_SOFTCLIP;
    MapOpts.flag |= MM_F_LONG_CIGAR;
    MapOpts.flag |= MM_F_EQX;
    MapOpts.flag |= MM_F_NO_PRINT_2ND;
    if (settings.NoTrimming) {
        MapOpts.flag |= MM_F_HARD_MLEVEL;
        MapOpts.mask_level = 0;
    }
    MapOpts.min_join_flank_ratio = 0.5;

    switch (settings.AlignMode) {
        case AlignmentMode::SUBREADS:
            *preset = "SUBREADS";
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
            *preset = "CCS";
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
            *preset = "ISOSEQ";
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
        case AlignmentMode::UNROLLED:
            *preset = "UNROLLED";
            MapOpts.flag |= MM_F_SPLICE | MM_F_SPLICE_FOR;
            MapOpts.max_gap = 10000;
            MapOpts.max_gap_ref = 2000;
            MapOpts.bw = 2000;
            MapOpts.a = 1;
            MapOpts.b = 2;
            MapOpts.q = 2;
            MapOpts.e = 1;
            MapOpts.q2 = 32;
            MapOpts.e2 = 0;
            MapOpts.zdrop = 200;
            MapOpts.zdrop_inv = 100;
            MapOpts.min_mid_occ = 100;
            MapOpts.min_dp_max = 200;
            MapOpts.noncan = 0;
            break;
        default:
            PBLOG_FATAL << "No AlignmentMode --preset selected!";
            std::exit(EXIT_FAILURE);
    }
    if (settings.GapOpen1 >= 0) MapOpts.q = settings.GapOpen1;
    if (settings.GapOpen2 >= 0) MapOpts.q2 = settings.GapOpen2;
    if (settings.GapExtension1 >= 0) MapOpts.e = settings.GapExtension1;
    if (settings.GapExtension2 >= 0) MapOpts.e2 = settings.GapExtension2;
    if (settings.MatchScore >= 0) MapOpts.a = settings.MatchScore;
    if (settings.MismatchPenalty >= 0) MapOpts.b = settings.MismatchPenalty;
    if (settings.Zdrop >= 0) MapOpts.zdrop = settings.Zdrop;
    if (settings.ZdropInv >= 0) MapOpts.zdrop_inv = settings.ZdropInv;
    if (settings.NonCanon >= 0) MapOpts.noncan = settings.NonCanon;
    if (settings.MaxIntronLength >= 0) mm_mapopt_max_intron_len(&MapOpts, settings.MaxIntronLength);
    if (settings.Bandwidth >= 0) MapOpts.bw = settings.Bandwidth;
    if (settings.NoSpliceFlank) MapOpts.flag &= ~MM_F_SPLICE_FLANK;
    if (settings.LongJoinFlankRatio >= 0)
        MapOpts.min_join_flank_ratio = settings.LongJoinFlankRatio;

    if ((MapOpts.q != MapOpts.q2 || MapOpts.e != MapOpts.e2) &&
        !(MapOpts.e > MapOpts.e2 && MapOpts.q + MapOpts.e < MapOpts.q2 + MapOpts.e2)) {
        PBLOG_FATAL << "Violation of dual gap penalties, E1>E2 and O1+E1<O2+E2";
        std::exit(EXIT_FAILURE);
    }

    if ((MapOpts.q + MapOpts.e) + (MapOpts.q2 + MapOpts.e2) > 127) {
        PBLOG_FATAL << "Violation of scoring system ({-O}+{-E})+({-O2}+{-E2}) <= 127";
        std::exit(EXIT_FAILURE);
    }

    if (MapOpts.zdrop < MapOpts.zdrop_inv) {
        PBLOG_FATAL << "Z-drop should not be less than inversion-Z-drop";
        std::exit(EXIT_FAILURE);
    }
}

void MM2Helper::PostInit(const MM2Settings& settings, const std::string& preset,
                         const bool postAlignParameter)
{
    mm_mapopt_update(&MapOpts, Idx->idx_);

    if (Idx->idx_->k <= 0 || Idx->idx_->w <= 0) {
        PBLOG_FATAL << "Index parameter -k and -w must be positive.";
        std::exit(EXIT_FAILURE);
    }

    PBLOG_DEBUG << "Minimap2 parameters based on preset: " << preset;
    PBLOG_DEBUG << "Kmer size              : " << Idx->idx_->k;
    PBLOG_DEBUG << "Minimizer window size  : " << Idx->idx_->w;
    PBLOG_DEBUG << "Homopolymer compressed : " << std::boolalpha
                << static_cast<bool>(Idx->idx_->flag & MM_I_HPC);
    if (postAlignParameter) {
        PBLOG_DEBUG << "Gap open 1             : " << MapOpts.q;
        PBLOG_DEBUG << "Gap open 2             : " << MapOpts.q2;
        PBLOG_DEBUG << "Gap extension 1        : " << MapOpts.e;
        PBLOG_DEBUG << "Gap extension 2        : " << MapOpts.e2;
        PBLOG_DEBUG << "Match score            : " << MapOpts.a;
        PBLOG_DEBUG << "Mismatch penalty       : " << MapOpts.b;
        PBLOG_DEBUG << "Z-drop                 : " << MapOpts.zdrop;
        PBLOG_DEBUG << "Z-drop inv             : " << MapOpts.zdrop_inv;
        PBLOG_DEBUG << "Bandwidth              : " << MapOpts.bw;
        PBLOG_DEBUG << "Long join flank ratio  : " << MapOpts.min_join_flank_ratio;
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
    auto tbuf = std::make_unique<ThreadBuffer>();
    auto result = std::make_unique<std::vector<AlignedRecord>>();
    result->reserve(records->size());

    for (auto& record : *records) {
        std::vector<AlignedRecord> localResults = Align(record, filter, tbuf);
        for (const auto& aln : localResults) {
            if (aln.IsAligned) {
                *alignedReads += 1;
                break;
            }
        }

        for (auto&& a : localResults)
            result->emplace_back(std::move(a));
    }

    return result;
}

std::vector<AlignedRecord> MM2Helper::Align(const BAM::BamRecord& record, const FilterFunc& filter,
                                            std::unique_ptr<ThreadBuffer>& tbuf) const
{
    std::vector<AlignedRecord> localResults;
    if (record.IsMapped() && record.Impl().IsSupplementaryAlignment()) return localResults;

    std::unique_ptr<ThreadBuffer> tbufLocal;
    if (!tbuf) tbufLocal = std::make_unique<ThreadBuffer>();

    int numAlns;
    const auto seq = record.Sequence(BAM::Orientation::NATIVE);
    std::unique_ptr<BAM::BamRecord> unalignedCopy;
    if (record.IsMapped()) {
        unalignedCopy = std::make_unique<BAM::BamRecord>(record.Header());
        unalignedCopy->Impl().SetSequenceAndQualities(
            seq, record.Qualities(BAM::Orientation::NATIVE).Fastq());
        unalignedCopy->Impl().SetReverseStrand(false);
        unalignedCopy->Impl().SetMapped(false);
        unalignedCopy->Impl().CigarData("");
        for (const auto& tag : record.Impl().Tags())
            if (tag.first != "RG") unalignedCopy->Impl().AddTag(tag.first, tag.second);
        unalignedCopy->ReadGroupId(record.ReadGroupId());
        unalignedCopy->Impl().Name(record.FullName());
    }

    const int qlen = seq.length();
    auto alns = mm_map(Idx->idx_, qlen, seq.c_str(), &numAlns,
                       tbufLocal ? tbufLocal->tbuf_ : tbuf->tbuf_, &MapOpts, nullptr);

    std::vector<int> used;
    std::vector<int32_t> queryHits(seq.size(), 0);

    for (int i = 0; i < numAlns; ++i) {
        auto& aln = alns[i];
        // if no alignment, continue
        if (aln.p == nullptr) continue;
        // secondary alignment
        if (aln.id != aln.parent) continue;
        used.emplace_back(i);
        if (alnMode_ == AlignmentMode::UNROLLED) break;
    }

    const auto SetQryHits = [&](mm_reg1_t& aln, int* begin, int* end) {
        bool started = false;
        bool ended = false;

        int l = aln.qs;
        int r = aln.qe;
        for (int s = l; s < r; ++s) {
            if (!started) {
                if (queryHits[s] == 0) {
                    queryHits[s] = 1;
                    *begin = s;
                    started = true;
                }
            } else {
                if (!ended) {
                    if (queryHits[s] == 0) {
                        queryHits[s] = 1;
                    } else if (queryHits[s] == 1) {
                        *end = s;
                        ended = true;
                        break;
                    }
                }
            }
        }
        if (!ended) *end = r;
        return started;
    };

    const auto AlignAndTrim = [&](const int idx, const bool trim) {
        if (maxNumAlns_ > 0 && static_cast<int32_t>(localResults.size()) >= maxNumAlns_) return;
        auto& aln = alns[idx];
        int begin = trim ? 0 : aln.qs;
        int end = trim ? 0 : aln.qe;
        if (trim && !SetQryHits(aln, &begin, &end)) return;
        const int32_t refId = aln.rid;
        const Data::Strand strand = aln.rev ? Data::Strand::REVERSE : Data::Strand::FORWARD;
        int refStartOffset = 0;
        Data::Cigar cigar;
        if (trim)
            cigar = RenderCigar(&aln, qlen, MapOpts.flag, begin, end, &refStartOffset);
        else
            cigar = RenderCigar(&aln, qlen, MapOpts.flag);
        const Data::Position refStart = aln.rs + refStartOffset;
        auto mapped = BAM::BamRecord::Mapped(unalignedCopy ? *unalignedCopy : record, refId,
                                             refStart, strand, cigar, aln.mapq);
        if (mapped.Impl().HasTag("rm")) mapped.Impl().RemoveTag("rm");
        mapped.Impl().SetSupplementaryAlignment(aln.sam_pri == 0);
        AlignedRecord alnRec{std::move(mapped)};
        if (filter(alnRec)) localResults.emplace_back(std::move(alnRec));
    };

    if (!trimRepeatedMatches_ || used.size() <= 1) {
        for (const auto i : used) {
            AlignAndTrim(i, false);
        }
    } else {
        std::vector<std::pair<int, int>> queryIntervals;
        std::vector<std::pair<int, int>> qryIdxOverlaps;
        std::vector<std::pair<int, int>> refIntervals;
        std::vector<std::pair<int, int>> refIdxOverlaps;

        for (const auto i : used) {
            auto& aln = alns[i];
            int l = aln.rev ? qlen - aln.qe : aln.qs;
            int r = aln.rev ? qlen - aln.qs : aln.qe;
            for (size_t j = 0; j < queryIntervals.size(); ++j) {
                const auto& s_e = queryIntervals[j];
                if (l <= s_e.second && r >= s_e.first) qryIdxOverlaps.emplace_back(i, used[j]);
            }
            queryIntervals.emplace_back(l, r);
            for (size_t j = 0; j < refIntervals.size(); ++j) {
                const auto& s_e = refIntervals[j];
                if (aln.rs <= s_e.second && aln.re >= s_e.first) {
                    if (aln.rs < alns[used[j]].rs)
                        refIdxOverlaps.emplace_back(i, used[j]);
                    else
                        refIdxOverlaps.emplace_back(used[j], i);
                }
            }
            refIntervals.emplace_back(aln.rs, aln.re);
        }
        std::vector<std::pair<int, int>> primaryAlignments;
        if (!refIdxOverlaps.empty()) {
            for (const auto& refIdx : refIdxOverlaps) {
                for (const auto& qryIdx : qryIdxOverlaps) {
                    if (refIdx.first == qryIdx.first || refIdx.first == qryIdx.second ||
                        refIdx.second == qryIdx.first || refIdx.second == qryIdx.second) {
                        primaryAlignments.emplace_back(refIdx);
                        break;
                    }
                }
            }
        }
        if (!qryIdxOverlaps.empty() && !primaryAlignments.empty()) {
            AlignAndTrim(primaryAlignments[0].first, true);
            AlignAndTrim(primaryAlignments[0].second, true);
            for (const auto i : used) {
                if (i == primaryAlignments[0].first) continue;
                if (i == primaryAlignments[0].second) continue;
                AlignAndTrim(i, true);
            }
            for (auto& l : localResults)
                l.Record.Impl().AddTag("rm", 1);
        } else {
            for (const auto i : used) {
                AlignAndTrim(i, true);
            }
        }
    }
    const auto numAlignments = localResults.size();
    if (numAlignments > 1) {
        std::vector<std::string> sas;
        std::vector<std::pair<int32_t, int32_t>> spans;
        for (size_t j = 0; j < numAlignments; ++j) {
            std::ostringstream sa;
            const auto& rec = localResults.at(j).Record;
            int32_t qryStart = BAM::IsCcsOrTranscript(rec.Type()) ? 0 : rec.QueryStart();
            const auto qqs = rec.AlignedStart() - qryStart;
            const auto qqe = rec.AlignedEnd() - qryStart;
            const auto qrs = rec.ReferenceStart();
            const auto qre = rec.ReferenceEnd();
            const bool qrev = rec.AlignedStrand() == Data::Strand::REVERSE;
            int l_M, l_I = 0, l_D = 0, clip5 = 0, clip3 = 0;
            if (qqe - qqs < qre - qrs)
                l_M = qqe - qqs, l_D = (qre - qrs) - l_M;
            else
                l_M = qre - qrs, l_I = (qqe - qqs) - l_M;
            clip5 = qrev ? qlen - qqe : qqs;
            clip3 = qrev ? qqs : qlen - qqe;
            sa << Idx->idx_->seq[rec.ReferenceId()].name << ',' << qrs + 1 << ',' << "+-"[qrev]
               << ',';
            if (clip5) sa << clip5 << 'S';
            if (l_M) sa << l_M << 'M';
            if (l_I) sa << l_I << 'I';
            if (l_D) sa << l_D << 'D';
            if (clip3) sa << clip3 << 'S';
            sa << ',' << static_cast<int>(rec.MapQuality()) << ',' << rec.NumMismatches() << ';';
            sas.emplace_back(sa.str());
            if (qrev)
                spans.emplace_back(clip3, qlen - clip5);
            else
                spans.emplace_back(clip5, qlen - clip3);
        }
        if (trimRepeatedMatches_)
            for (size_t i = 0; i < numAlignments; ++i) {
                int32_t fixedBegin = spans[i].first;
                int32_t fixedEnd = spans[i].second;
                for (size_t j = i + 1; j < numAlignments; ++j) {
                    int32_t curBegin = spans[j].first;
                    int32_t curEnd = spans[j].second;
                    if (fixedEnd > curBegin && fixedBegin < curEnd) {
                        PBLOG_FATAL << "Overlapping intervals: " << record.FullName() << "\t"
                                    << fixedBegin << "-" << fixedEnd << "\t" << curBegin << "-"
                                    << curEnd;
                    }
                }
            }
        for (size_t i = 0; i < numAlignments; ++i) {
            std::ostringstream sa;
            for (size_t j = 0; j < numAlignments; ++j) {
                if (i == j) continue;
                sa << sas[j];
            }
            const auto sastr = sa.str();
            if (!sastr.empty()) {
                if (localResults[i].Record.Impl().HasTag("SA"))
                    localResults[i].Record.Impl().EditTag("SA", sastr);
                else
                    localResults[i].Record.Impl().AddTag("SA", sastr);
            }
        }
    }

    // cleanup
    for (int i = 0; i < numAlns; ++i)
        if (alns[i].p) free(alns[i].p);
    free(alns);

    if (localResults.empty()) {
        if (record.IsMapped()) {
            const auto RemovePbmm2MappedTags = [](BAM::BamRecord& r) {
                for (const auto& t : {"SA", "rm", "mc"})
                    if (r.Impl().HasTag(t)) r.Impl().RemoveTag(t);
            };
            if (unalignedCopy) {
                RemovePbmm2MappedTags(*unalignedCopy);
                localResults.emplace_back(*unalignedCopy);
            }
        } else {
            localResults.emplace_back(AlignedRecord{record});
        }
    }

    return localResults;
}

std::vector<AlignedRecord> MM2Helper::Align(const BAM::BamRecord& record) const
{
    auto tbuf = std::make_unique<ThreadBuffer>();
    const auto noopFilter = [](const AlignedRecord&) { return true; };
    return Align(record, noopFilter, tbuf);
}

std::vector<AlignedRecord> MM2Helper::Align(const BAM::BamRecord& record,
                                            const FilterFunc& filter) const
{
    auto tbuf = std::make_unique<ThreadBuffer>();
    return Align(record, filter, tbuf);
}

std::vector<AlignedRecord> MM2Helper::Align(const BAM::BamRecord& record,
                                            std::unique_ptr<ThreadBuffer>& tbuf) const
{
    const auto noopFilter = [](const AlignedRecord&) { return true; };
    return Align(record, noopFilter, tbuf);
}

std::vector<BAM::SequenceInfo> MM2Helper::SequenceInfos() const { return Idx->SequenceInfos(); }

Index::Index(const std::vector<BAM::FastaSequence>& refs, const mm_idxopt_t& opts)
{
    IndexFrom(refs, opts);
}

Index::Index(std::vector<BAM::FastaSequence>&& refs, const mm_idxopt_t& opts)
    : refs_{std::move(refs)}
{
    IndexFrom(refs_, opts);
}

void Index::IndexFrom(const std::vector<BAM::FastaSequence>& refs, const mm_idxopt_t& opts)
{
    const auto numRefs = refs.size();
    seq_ = (const char**)calloc(numRefs + 1, sizeof(char*));
    for (size_t i = 0; i < numRefs; ++i)
        seq_[i] = refs[i].Bases().c_str();
    name_ = (const char**)calloc(numRefs + 1, sizeof(char*));
    for (size_t i = 0; i < numRefs; ++i)
        name_[i] = refs[i].Name().c_str();

    idx_ = mm_idx_str(opts.w, opts.k, opts.flag & MM_I_HPC, 0, numRefs, seq_, name_);
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

Index::~Index()
{
    free(seq_);
    free(name_);

    mm_idx_destroy(idx_);
}

std::vector<BAM::SequenceInfo> Index::SequenceInfos() const
{
    std::vector<BAM::SequenceInfo> result;
    for (unsigned i = 0; i < idx_->n_seq; ++i) {
        const std::string name = idx_->seq[i].name;
        const std::string len = std::to_string(idx_->seq[i].len);
        result.emplace_back(BAM::SequenceInfo(name, len));
    }
    return result;
}

AlignedRecord::AlignedRecord(BAM::BamRecord record) : Record(std::move(record))
{
    IsAligned = Record.IsMapped();
    if (IsAligned) ComputeAccuracyBases();
}

void AlignedRecord::ComputeAccuracyBases()
{
    int32_t ins = 0;
    int32_t del = 0;
    int32_t mismatch = 0;
    int32_t match = 0;
    for (const auto& cigar : Record.CigarData()) {
        int32_t len = cigar.Length();
        switch (cigar.Type()) {
            case Data::CigarOperationType::INSERTION:
                ins += len;
                break;
            case Data::CigarOperationType::DELETION:
                del += len;
                break;
            case Data::CigarOperationType::SEQUENCE_MISMATCH:
                mismatch += len;
                break;
            case Data::CigarOperationType::REFERENCE_SKIP:
                break;
            case Data::CigarOperationType::SEQUENCE_MATCH:
            case Data::CigarOperationType::ALIGNMENT_MATCH:
                match += len;
                break;
            case Data::CigarOperationType::PADDING:
            case Data::CigarOperationType::SOFT_CLIP:
            case Data::CigarOperationType::HARD_CLIP:
                break;
            case Data::CigarOperationType::UNKNOWN_OP:
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
    if (Record.Impl().HasTag("mc"))
        Record.Impl().EditTag("mc", static_cast<float>(Concordance));
    else
        Record.Impl().AddTag("mc", static_cast<float>(Concordance));
}

}  // namespace minimap2
}  // namespace PacBio
