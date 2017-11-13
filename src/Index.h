#pragma once

#include <map>
#include <vector>

#include <pbbam/FastaReader.h>
#include <pbbam/MD5.h>
#include <pbbam/SequenceInfo.h>

#include <minimap.h>

extern "C" void mm_idxopt_init(mm_idxopt_t*);

namespace PacBio {
namespace minimap2 {

struct IndexOptions
{
    IndexOptions(int nthreads = 3) : nthreads_{nthreads} {
        // not included in minimap.h for some reason
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
        PacBio::BAM::FastaReader faRdr(fname);
        PacBio::BAM::FastaSequence rec;
        while (faRdr.GetNext(rec))
            m5m_[rec.Name()] = PacBio::BAM::MD5Hash(rec.Bases());
    }

    ~Index() {
        mm_idx_destroy(idx_);
    }

    std::vector<PacBio::BAM::SequenceInfo> SequenceInfos() const {
        using PacBio::BAM::SequenceInfo;
        std::vector<SequenceInfo> result;
        for (unsigned i = 0; i < idx_->n_seq; ++i) {
            const std::string name = idx_->seq[i].name;
            const std::string len = std::to_string(idx_->seq[i].len);
            const std::string md5 = m5m_.at(name);
            result.emplace_back(SequenceInfo(name, len).Checksum(md5));
        }
        return result;
    }

    mm_idx_t* idx_;
    std::map<std::string, std::string> m5m_;
};

}}
