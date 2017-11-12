
#include <iostream>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>

#include <minimap.h>

using namespace PacBio::BAM;

// not included in minimap.h for some reason
extern "C" void mm_idxopt_init(mm_idxopt_t*);

struct IndexOptions
{
    IndexOptions(int nthreads = 3) : nthreads_{nthreads} {
        mm_idxopt_init(&opts_);
        // set default map-pb options
        // map-pb settings (aka -Hk19)
		opts_.is_hpc = 1;  // homopolymer compression (-H)
        opts_.k = 19;  // k-mer of 19 (-k19)
		opts_.batch_size = 0x7fffffffffffffffL;  // always build a uni-part index
    }

    mm_idxopt_t opts_;
    int nthreads_;
};

struct Index
{
    Index(const std::string& fname, const IndexOptions& opts) : idx_{nullptr} {
        auto rdr = mm_idx_reader_open(fname.c_str(), &opts.opts_, nullptr);
        if (!rdr) throw std::runtime_error("unable to load reference for indexing!");
        idx_ = mm_idx_reader_read(rdr, opts.nthreads_);
        if (!idx_) throw std::runtime_error("unable to index reference!");
        mm_idx_reader_close(rdr);
    }

    ~Index() {
        mm_idx_destroy(idx_);
    }

    std::vector<SequenceInfo> SequenceInfos() const {
        std::vector<SequenceInfo> result;
        for (unsigned i = 0; i < idx_->n_seq; ++i)
            result.emplace_back(SequenceInfo(std::string(idx_->seq[i].name),
                                             std::to_string(idx_->seq[i].len)));
        return result;
    }

    mm_idx_t* idx_;
};

struct MapOptions
{
    MapOptions() {
		mm_mapopt_init(&opts_);
		opts_.flag |= 4;  // always perform alignment
    }

    void Update(const Index& idx) {
        mm_mapopt_update(&opts_, idx.idx_);
    }

    mm_mapopt_t opts_;
};

struct ThreadBuffer
{
    ThreadBuffer() {
        tbuf_ = mm_tbuf_init();
    }

    ~ThreadBuffer() {
        mm_tbuf_destroy(tbuf_);
    }

    mm_tbuf_t* tbuf_;
};

Cigar RenderCigar(const mm_reg1_t* const r, const int qlen, const int opt_flag)
{
    Cigar cigar;

    if (r->p == nullptr)
        return cigar;

    uint32_t k, clip_len[2];
    clip_len[0] = r->rev? qlen - r->qe : r->qs;
    clip_len[1] = r->rev? r->qs : qlen - r->qe;
    const char clip_char = 'S'; /* (sam_flag & 0x800) && !(opt_flag & MM_F_SOFTCLIP) ? 'H' : 'S'; */

    if (clip_len[0]) cigar.emplace_back(clip_char, clip_len[0]);
    for (k = 0; k < r->p->n_cigar; ++k)
        cigar.emplace_back("MIDN=X"[r->p->cigar[k] & 0xf], r->p->cigar[k] >> 4);
    if (clip_len[1]) cigar.emplace_back(clip_char, clip_len[1]);

#if 0
    k = 0;
    for (const auto op : cigar)
        if (op.Type() != CigarOperationType::DELETION && op.Type() != CigarOperationType::REFERENCE_SKIP)
            k += op.Length();

    if (static_cast<int>(k) != qlen)
        std::cerr << "cigar: " << k << ", qlen: " << qlen << std::endl;

    std::cerr << "cigar: " << cigar.ToStdString() << std::endl;
#endif

    return cigar;
}

std::vector<BamRecord> Align(const BamRecord& rec, const Index& idx, const MapOptions& mapOpts)
{
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
            std::cerr << "alignment for " << rec.FullName() << " exceeds maximum cigar operations (65535)" << std::endl;
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

int main(const int argc, const char** argv)
{
    IndexOptions idxOpts;
    MapOptions mapOpts;

    if (argc != 4) {
        std::cerr << "usage: pbmm2 REFFA QUERYBAM ALNOUT" << std::endl;
        return -1;
    }

    std::string refFile(argv[1]);
    std::string qryFile(argv[2]);
    std::string alnFile(argv[3]);

    Index idx(refFile, idxOpts);
    mapOpts.Update(idx);

    BamReader qryRdr(qryFile);
    BamHeader hdr(qryRdr.Header());

    {
        for (const auto si : idx.SequenceInfos())
            hdr.AddSequence(si);
        hdr.AddProgram(ProgramInfo("pbmm2").Name("pbmm2").Version("0.0.1"));
    }

    BamWriter out(alnFile, hdr);

    BamRecord rec;
    while (qryRdr.GetNext(rec)) {
        for (const auto& aln : Align(rec, idx, mapOpts))
            out.Write(aln);
    }

    return 0;
}
