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

static void SetGlobalLogger(const Logging::LogLevel logLevel = Logging::LogLevel::FATAL)
{
    PacBio::Logging::InstallSignalHandlers(
        PacBio::Logging::Logger::Default(new PacBio::Logging::Logger(std::cerr, logLevel)));
}

static std::vector<BAM::BamRecord> SimpleAlign()
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

TEST(MM2Test, SimpleAlign)
{
    SetGlobalLogger();
    std::vector<BAM::BamRecord> alignedBam = SimpleAlign();
    EXPECT_EQ(96, alignedBam.size());
}

TEST(MM2Test, Filter)
{
    SetGlobalLogger();
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
    EXPECT_EQ(85, alignedBam.size());
}

TEST(MM2Test, UseCommonThreadBuffer)
{
    SetGlobalLogger();
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
    EXPECT_EQ(96, alignedBam.size());
}

TEST(MM2Test, FilterAndBuffer)
{
    SetGlobalLogger();
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
    EXPECT_EQ(85, alignedBam.size());
}

TEST(MM2Test, AlignCCS)
{
    SetGlobalLogger();
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
    EXPECT_EQ(8, alignedBam.size());
    EXPECT_EQ(11501, alignedBases);
}

TEST(MM2Test, AlignCCSWithOverrides)
{
    SetGlobalLogger(Logging::LogLevel::DEBUG);
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
    EXPECT_EQ(9, alignedBam.size());
    EXPECT_EQ(11704, alignedBases);
}

// For a proper test of idempotence, see tests/cram/idempotence.t
TEST(MM2Test, ReAlign)
{
    SetGlobalLogger();
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

    EXPECT_EQ(96, alignedBam.size());
    EXPECT_EQ(96, idempotent.size());
}

static std::vector<BAM::BamRecord> FastaRefAlign()
{
    const auto refFile = tests::DataDir + '/' + "ecoliK12_pbi_March2013.fasta";

    MM2Settings settings;
    MM2Helper mm2helper(BAM::FastaReader::ReadAll(refFile), settings);
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
    SetGlobalLogger();

    std::vector<BAM::BamRecord> alignedFastaRefBam = FastaRefAlign();
    EXPECT_EQ(96, alignedFastaRefBam.size());

    std::vector<BAM::BamRecord> alignedRefFileBam = SimpleAlign();
    EXPECT_EQ(96, alignedRefFileBam.size());

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

}  // namespace MM2Tests
}  // namespace PacBio
