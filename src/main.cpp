#include <iostream>
#include <memory>
#include <thread>
#include <vector>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>

#include "Index.h"
#include "Mapping.h"
#include "WorkQueue.h"

using namespace PacBio::BAM;
using namespace PacBio::Parallel;
using namespace PacBio::minimap2;

typedef std::vector<BamRecord> Results;

void WriteRecords(BamWriter& out, Results&& results)
{
    for (const auto& aln : results)
        out.Write(aln);
}

void WriterThread(WorkQueue<Results>& queue, std::unique_ptr<BamWriter> out)
{
    while (queue.ConsumeWith(WriteRecords, std::ref(*out)));
}

int main(const int argc, const char** argv)
{
    using std::cref;
    using std::move;
    using std::ref;

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

    const int numThreads = std::thread::hardware_concurrency();
    WorkQueue<Results> queue(numThreads);

    std::unique_ptr<BamWriter> out(new BamWriter(alnFile, hdr));
    std::future<void> writer = std::async(std::launch::async, WriterThread, ref(queue), move(out));

    BamRecord rec;
    while (qryRdr.GetNext(rec)) queue.ProduceWith(&Align, rec, cref(idx), cref(mapOpts));

    queue.Finalize();
    writer.wait();

    return 0;
}
