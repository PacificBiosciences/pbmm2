// Author: Armin Töpfer

#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaReader.h>
#include <pbcopper/logging/Logging.h>

#include <pbmm2/AlignmentMode.h>
#include <pbmm2/MM2Helper.h>

#include <TestData.h>

namespace PacBio {
namespace MM2Tests {

using namespace PacBio::minimap2;

static std::vector<BAM::BamRecord> SimpleAlignBAM()
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<BAM::BamRecord> alignedBam;
    for (const auto& record : reader) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record);
        for (const auto& aln : alignments)
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
    }
    return alignedBam;
}

static std::vector<CompatMappedRead> SimpleAlignRead()
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<CompatMappedRead> results;
    for (const auto& record : reader) {
        const auto inputRead = record.ToRead();
        const std::vector<AlignedRead> alignments = mm2helper.Align(inputRead);
        for (const auto& aln : alignments)
            if (aln.IsAligned) results.emplace_back(std::move(aln.Record));
    }
    return results;
}

TEST(MM2Test, SimpleAlignBAM)
{
    std::vector<BAM::BamRecord> alignedBam = SimpleAlignBAM();
    EXPECT_EQ(96ul, alignedBam.size());
}

TEST(MM2Test, SimpleAlignRead)
{
    std::vector<CompatMappedRead> alignedReads = SimpleAlignRead();
    EXPECT_EQ(96ul, alignedReads.size());
}

TEST(MM2Test, VectorAlignRead)
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    auto rawReads = std::make_unique<std::vector<Data::Read>>();
    int32_t numReadsBam = 0;
    for (const auto& record : reader) {
        ++numReadsBam;
        rawReads->emplace_back(record.ToRead());
    }

    int32_t alignedReads{};
    const auto myFilter = [](const AlignedRead& r) { return true; };
    const auto alns = mm2helper.Align(rawReads, myFilter, &alignedReads);
    EXPECT_EQ(numReadsBam, alignedReads);
}

TEST(MM2Test, FilterBAM)
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<BAM::BamRecord> alignedBam;

    const auto myFilter = [](const AlignedRecord& r) { return r.Concordance > 90; };

    for (const auto& record : reader) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record, myFilter);
        for (const auto& aln : alignments) {
            EXPECT_TRUE(aln.Concordance > 90.0);
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
        }
    }
    EXPECT_EQ(85ul, alignedBam.size());
}

TEST(MM2Test, FilterRead)
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<CompatMappedRead> alignedReads;

    const auto myFilter = [](const AlignedRead& r) { return r.Concordance > 90; };

    for (const auto& record : reader) {
        const auto inputRead = record.ToRead();
        const std::vector<AlignedRead> alignments = mm2helper.Align(inputRead, myFilter);
        for (const auto& aln : alignments) {
            EXPECT_TRUE(aln.Concordance > 90.0);
            if (aln.IsAligned) alignedReads.emplace_back(std::move(aln.Record));
        }
    }
    EXPECT_EQ(85ul, alignedReads.size());
}

TEST(MM2Test, UseCommonThreadBuffer)
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<BAM::BamRecord> alignedBam;

    // Use _one_ thread buffer per thread, do not reuse between threads!
    auto tbuf = std::make_unique<ThreadBuffer>();
    for (const auto& record : reader) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record, tbuf);
        for (const auto& aln : alignments) {
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
        }
    }
    EXPECT_EQ(96ul, alignedBam.size());
}

TEST(MM2Test, FilterAndBuffer)
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<BAM::BamRecord> alignedBam;

    const auto myFilter = [](const AlignedRecord& r) { return r.Concordance > 90; };
    auto tbuf = std::make_unique<ThreadBuffer>();

    for (const auto& record : reader) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record, myFilter, tbuf);
        for (const auto& aln : alignments) {
            EXPECT_TRUE(aln.Concordance > 90.0);
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
        }
    }
    EXPECT_EQ(85ul, alignedBam.size());
}

TEST(MM2Test, AlignCCS)
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    settings.AlignMode = AlignmentMode::CCS;

    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "m54075_180905_221350.ccs.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<BAM::BamRecord> alignedBam;
    int32_t alignedBases = 0;
    for (const auto& record : reader) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record);
        for (const auto& aln : alignments) {
            alignedBases += aln.NumAlignedBases;
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
        }
    }
    EXPECT_EQ(8ul, alignedBam.size());
    EXPECT_EQ(11501, alignedBases);
}

TEST(MM2Test, AlignCCSWithOverrides)
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    settings.AlignMode = AlignmentMode::CCS;

    settings.Kmer = 10;
    settings.MinimizerWindowSize = 5;

    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "m54075_180905_221350.ccs.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<BAM::BamRecord> alignedBam;
    int32_t alignedBases = 0;
    for (const auto& record : reader) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record);
        for (const auto& aln : alignments) {
            alignedBases += aln.NumAlignedBases;
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
        }
    }
    EXPECT_EQ(9ul, alignedBam.size());
    EXPECT_EQ(11704, alignedBases);
}

// For a proper test of idempotence, see tests/cram/idempotence.t
TEST(MM2Test, ReAlignBAM)
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<BAM::BamRecord> alignedBam;
    for (const auto& record : reader) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record);
        for (const auto& aln : alignments)
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
    }

    std::vector<BAM::BamRecord> idempotent;
    for (const auto& record : alignedBam) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record);
        for (const auto& aln : alignments)
            idempotent.emplace_back(std::move(aln.Record));
    }

    EXPECT_EQ(96ul, alignedBam.size());
    EXPECT_EQ(96ul, idempotent.size());
}

TEST(MM2Test, ReAlignRead)
{
    // Data::Read and friends never flip the Sequence or PulseWidths

    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";
    MM2Settings settings;
    MM2Helper mm2helper(refFile, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<CompatMappedRead> alignedReads;
    for (const auto& record : reader) {
        const auto inputRead = record.ToRead();
        const std::vector<AlignedRead> alignments = mm2helper.Align(inputRead);
        for (const auto& aln : alignments)
            if (aln.IsAligned) alignedReads.emplace_back(std::move(aln.Record));
    }

    std::vector<CompatMappedRead> idempotent;
    for (const auto& record : alignedReads) {
        // need to skip aligned supplementary alignments
        if (record.IsSupplementaryAlignment()) continue;

        const std::vector<AlignedRead> alignments = mm2helper.Align(record);
        for (const auto& aln : alignments)
            idempotent.emplace_back(std::move(aln.Record));
    }

    EXPECT_EQ(96ul, alignedReads.size());
    EXPECT_EQ(96ul, idempotent.size());
}

static std::vector<BAM::BamRecord> FastaRefAlign()
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";

    MM2Settings settings;
    const std::vector<BAM::FastaSequence> refs = BAM::FastaReader::ReadAll(refFile);
    MM2Helper mm2helper(refs, settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<BAM::BamRecord> alignedBam;
    for (const auto& record : reader) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record);
        for (const auto& aln : alignments)
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
    }
    return alignedBam;
}

TEST(MM2Test, FastaRefAlign)
{

    std::vector<BAM::BamRecord> alignedFastaRefBam = FastaRefAlign();
    EXPECT_EQ(96ul, alignedFastaRefBam.size());

    std::vector<BAM::BamRecord> alignedRefFileBam = SimpleAlignBAM();
    EXPECT_EQ(96ul, alignedRefFileBam.size());

    for (size_t i = 0; i < alignedFastaRefBam.size(); ++i) {
        const auto& fromFasta = alignedFastaRefBam[i];
        const auto& fromFile = alignedRefFileBam[i];
        EXPECT_EQ(fromFasta.FullName(), fromFile.FullName());
        EXPECT_EQ(fromFasta.Sequence(), fromFile.Sequence());
        EXPECT_EQ(fromFasta.ReferenceId(), fromFile.ReferenceId());
        EXPECT_EQ(fromFasta.ReferenceStart(), fromFile.ReferenceStart());
        EXPECT_EQ(fromFasta.ReferenceEnd(), fromFile.ReferenceEnd());
        EXPECT_EQ(fromFasta.QueryStart(), fromFile.QueryStart());
        EXPECT_EQ(fromFasta.QueryEnd(), fromFile.QueryEnd());
        EXPECT_EQ(fromFasta.AlignedStart(), fromFile.AlignedStart());
        EXPECT_EQ(fromFasta.AlignedEnd(), fromFile.AlignedEnd());
        EXPECT_EQ(fromFasta.AlignedStrand(), fromFile.AlignedStrand());
        EXPECT_EQ(fromFasta.CigarData().ToStdString(), fromFile.CigarData().ToStdString());
        EXPECT_EQ(fromFasta.Impl().IsSupplementaryAlignment(),
                  fromFile.Impl().IsSupplementaryAlignment());
    }
}

static std::vector<BAM::BamRecord> FastaRefAlignMoveBAM()
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";

    MM2Settings settings;
    const std::vector<BAM::FastaSequence> refs = BAM::FastaReader::ReadAll(refFile);
    MM2Helper mm2helper(std::move(refs), settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<BAM::BamRecord> alignedBam;
    for (const auto& record : reader) {
        const std::vector<AlignedRecord> alignments = mm2helper.Align(record);
        for (const auto& aln : alignments)
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
    }
    return alignedBam;
}

TEST(MM2Test, FastaRefAlignMoveBAM)
{

    std::vector<BAM::BamRecord> alignedFastaRefBam = FastaRefAlignMoveBAM();
    EXPECT_EQ(96ul, alignedFastaRefBam.size());

    std::vector<BAM::BamRecord> alignedRefFileBam = SimpleAlignBAM();
    EXPECT_EQ(96ul, alignedRefFileBam.size());

    for (size_t i = 0; i < alignedFastaRefBam.size(); ++i) {
        const auto& fromFasta = alignedFastaRefBam[i];
        const auto& fromFile = alignedRefFileBam[i];
        EXPECT_EQ(fromFasta.FullName(), fromFile.FullName());
        EXPECT_EQ(fromFasta.Sequence(), fromFile.Sequence());
        EXPECT_EQ(fromFasta.ReferenceId(), fromFile.ReferenceId());
        EXPECT_EQ(fromFasta.ReferenceStart(), fromFile.ReferenceStart());
        EXPECT_EQ(fromFasta.ReferenceEnd(), fromFile.ReferenceEnd());
        EXPECT_EQ(fromFasta.QueryStart(), fromFile.QueryStart());
        EXPECT_EQ(fromFasta.QueryEnd(), fromFile.QueryEnd());
        EXPECT_EQ(fromFasta.AlignedStart(), fromFile.AlignedStart());
        EXPECT_EQ(fromFasta.AlignedEnd(), fromFile.AlignedEnd());
        EXPECT_EQ(fromFasta.AlignedStrand(), fromFile.AlignedStrand());
        EXPECT_EQ(fromFasta.CigarData().ToStdString(), fromFile.CigarData().ToStdString());
        EXPECT_EQ(fromFasta.Impl().IsSupplementaryAlignment(),
                  fromFile.Impl().IsSupplementaryAlignment());
    }
}

static std::vector<CompatMappedRead> FastaRefAlignMoveRead()
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";

    MM2Settings settings;
    const std::vector<BAM::FastaSequence> refs = BAM::FastaReader::ReadAll(refFile);
    MM2Helper mm2helper(std::move(refs), settings);
    const auto alnFile = tests::DataDir + '/' + "median.bam";
    BAM::EntireFileQuery reader(alnFile);
    std::vector<CompatMappedRead> alignedBam;
    for (const auto& record : reader) {
        const auto inputRead = record.ToRead();
        const std::vector<AlignedRead> alignments = mm2helper.Align(inputRead);
        for (const auto& aln : alignments)
            if (aln.IsAligned) alignedBam.emplace_back(std::move(aln.Record));
    }
    return alignedBam;
}

TEST(MM2Test, FastaRefAlignMoveRead)
{

    auto alignedFastaRefRead = FastaRefAlignMoveRead();
    EXPECT_EQ(96ul, alignedFastaRefRead.size());

    // text mixture of Data::Read and BAM
    auto alignedRefFileBam = SimpleAlignBAM();
    EXPECT_EQ(96ul, alignedRefFileBam.size());

    for (size_t i = 0; i < alignedFastaRefRead.size(); ++i) {
        const auto& fromFasta = alignedFastaRefRead[i];
        const auto& fromFile = alignedRefFileBam[i];
        EXPECT_EQ(fromFasta.FullName(), fromFile.FullName());
        EXPECT_EQ(fromFasta.Sequence(), fromFile.Sequence());
        EXPECT_EQ(fromFasta.ReferenceId(), fromFile.ReferenceId());
        EXPECT_EQ(fromFasta.ReferenceStart(), fromFile.ReferenceStart());
        EXPECT_EQ(fromFasta.ReferenceEnd(), fromFile.ReferenceEnd());
        EXPECT_EQ(fromFasta.QueryStart(), fromFile.QueryStart());
        EXPECT_EQ(fromFasta.QueryEnd(), fromFile.QueryEnd());
        EXPECT_EQ(fromFasta.AlignedStart(), fromFile.AlignedStart());
        EXPECT_EQ(fromFasta.AlignedEnd(), fromFile.AlignedEnd());
        EXPECT_EQ(fromFasta.AlignedStrand(), fromFile.AlignedStrand());
        EXPECT_EQ(fromFasta.CigarData().ToStdString(), fromFile.CigarData().ToStdString());
        EXPECT_EQ(fromFasta.Impl().IsSupplementaryAlignment(),
                  fromFile.Impl().IsSupplementaryAlignment());
    }
}

}  // namespace MM2Tests
}  // namespace PacBio
