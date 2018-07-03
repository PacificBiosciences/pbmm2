#pragma once

#include <map>
#include <vector>

#include <boost/algorithm/string.hpp>

#include <pbbam/FastaReader.h>
#include <pbbam/MD5.h>
#include <pbbam/SequenceInfo.h>

#include <minimap.h>

extern "C" void mm_idxopt_init(mm_idxopt_t*);

namespace PacBio {
namespace minimap2 {

struct IndexOptions
{
    IndexOptions(const int k = 19, const int w = 10, const bool hpc = true, const int nthreads = 3)
        : nthreads_{nthreads}
    {
        // not included in minimap.h for some reason
        mm_idxopt_init(&opts_);
        // set default map-pb options
        // map-pb settings (aka -Hk19)
        if (hpc) opts_.flag |= MM_I_HPC;         // homopolymer compression (-H)
        opts_.k = k;                             // k-mer of 19 (-k19)
        opts_.w = w;                             // minimizer window of 10 (default)
        opts_.batch_size = 0x7fffffffffffffffL;  // always build a uni-part index
    }

    mm_idxopt_t opts_;
    int nthreads_;
};

struct Index
{
    Index(const std::string& fname, const IndexOptions& opts) : idx_{nullptr}
    {
        using boost::algorithm::trim_right_if;
        using boost::algorithm::is_any_of;
        auto rdr = mm_idx_reader_open(fname.c_str(), &opts.opts_, nullptr);
        if (!rdr) throw std::runtime_error("unable to load reference for indexing!");
        idx_ = mm_idx_reader_read(rdr, opts.nthreads_);
        if (!idx_) throw std::runtime_error("unable to index reference!");
        mm_idx_reader_close(rdr);
        PacBio::BAM::FastaReader faRdr(fname);
        PacBio::BAM::FastaSequence rec;
        while (faRdr.GetNext(rec)) {
            std::string name = rec.Name();
            trim_right_if(name, is_any_of("\r\n"));
            std::string m5 = PacBio::BAM::MD5Hash(rec.Bases());
            m5m_[name] = m5;
            name = name.substr(0, name.find_first_of(" \f\n\r\t\v"));
            m5m_[name] = m5;
        }
    }

    ~Index() { mm_idx_destroy(idx_); }

    std::vector<PacBio::BAM::SequenceInfo> SequenceInfos() const
    {
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
}  // namespace minimap2
}  // namespace PacBio