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
    const auto myFilter = [](const AlignedRead&) { return true; };
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

void TestTrimCigarFlanks(std::string given, std::string expected, int expectedOffset)
{
    Data::Cigar c(given);
    int offset = 0;
    TrimCigarFlanks(&c, &offset);
    EXPECT_EQ(expected, c.ToStdString()) << " where cigar was " << given;
    EXPECT_EQ(offset, expectedOffset) << " where cigar was " << given;
}

TEST(MM2Test, TrimCigarFlanks)
{
    TestTrimCigarFlanks("3I4=", "3S4=", 0);
    TestTrimCigarFlanks("3D4=", "4=", 3);
    TestTrimCigarFlanks("3I4=5I", "3S4=5S", 0);
    TestTrimCigarFlanks("3D4=5D", "4=", 3);
    TestTrimCigarFlanks("2S3I4=", "5S4=", 0);
    TestTrimCigarFlanks("2S3D4=", "2S4=", 3);
    TestTrimCigarFlanks("2S3I4=5I6S", "5S4=11S", 0);
    TestTrimCigarFlanks("2S3D4=5D6S", "2S4=6S", 3);
    TestTrimCigarFlanks("2D3I", "3S", 2);
    TestTrimCigarFlanks("3I2D", "3S", 2);
    TestTrimCigarFlanks("1=2D3I", "1=3S", 0);
    TestTrimCigarFlanks("1=3I2D", "1=3S", 0);

    TestTrimCigarFlanks("2X4=", "2S4=", 2);
    TestTrimCigarFlanks("2X1D4=", "2S4=", 3);
    TestTrimCigarFlanks("2X1I4=", "3S4=", 2);
    TestTrimCigarFlanks("2X1D3I4=", "5S4=", 3);

    TestTrimCigarFlanks("3I2X4=", "5S4=", 2);
    TestTrimCigarFlanks("3I1D2X4=", "5S4=", 3);
    TestTrimCigarFlanks("3I2X1D4=", "5S4=", 3);

    TestTrimCigarFlanks("3I2X1D4=3X", "5S4=3S", 3);
    TestTrimCigarFlanks("3I2X1D4=3X2D", "5S4=3S", 3);
    TestTrimCigarFlanks("3I2X1D4=3X2D1I", "5S4=4S", 3);

    TestTrimCigarFlanks("3I2X1D4=3I", "5S4=3S", 3);
    TestTrimCigarFlanks("3I2X1D4=3I2X", "5S4=5S", 3);
    TestTrimCigarFlanks("3I2X1D4=3I2X1D", "5S4=5S", 3);

    TestTrimCigarFlanks("3I2X1D4=3D", "5S4=", 3);
    TestTrimCigarFlanks("3I2X1D4=3D1X", "5S4=1S", 3);
    TestTrimCigarFlanks("3I2X1D4=3D1X2I", "5S4=3S", 3);
    TestTrimCigarFlanks("3I2X1D4=3D1X2I1X", "5S4=4S", 3);

    TestTrimCigarFlanks("3I2X1D4D1X2I1X", "9S", 9);

    TestTrimCigarFlanks("10S1D1X1I1D1X1I10=1D1X1I1D1X1I10S", "14S10=14S", 4);

    // clang-format off
    TestTrimCigarFlanks(
        "913=1X437=1I159=1X105=1X308=1X57=2I131=1X25=2X84=1X46=20I215=1X248=1X136=1X152=1X21=1X9=1X13=1X32=436I8090S",
        "913=1X437=1I159=1X105=1X308=1X57=2I131=1X25=2X84=1X46=20I215=1X248=1X136=1X152=1X21=1X9=1X13=1X32=8526S",
        0);
    // clang-format on
}

}  // namespace MM2Tests
}  // namespace PacBio
