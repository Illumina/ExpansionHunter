//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "region_analysis/RegionAnalyzer.hh"

#include "gtest/gtest.h"

#include "input/GraphBlueprint.hh"
#include "input/RegionGraph.hh"
#include "region_spec/LocusSpecification.hh"

using namespace ehunter;

using graphtools::Graph;
using std::string;
using std::vector;

class AlignerTests : public ::testing::TestWithParam<std::string>
{
};

// TODO: Throw an error if there are no valid extensions?
TEST_P(AlignerTests, RegionAnalysis_ShortSingleUnitRepeat_Genotyped)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    vector<GenomicRegion> referenceRegions = { GenomicRegion(1, 1, 2) };

    NodeToRegionAssociation dummyAssociation;
    LocusSpecification regionSpec("region", referenceRegions, AlleleCount::kTwo, graph, dummyAssociation);
    VariantClassification classification(VariantType::kRepeat, VariantSubtype::kCommonRepeat);
    regionSpec.addVariantSpecification("repeat", classification, GenomicRegion(1, 1, 2), { 1 }, 1);

    SampleParameters sampleParams("dummy_sample", Sex::kFemale, 10, 5.0);
    HeuristicParameters heuristicParams(1000, 20, true, GetParam(), 4, 1, 5);

    graphtools::BlankAlignmentWriter blankAlignmentWriter;
    RegionAnalyzer regionAnalyzer(regionSpec, heuristicParams, blankAlignmentWriter);

    regionAnalyzer.processMates(
        Read(ReadId("read1", MateNumber::kFirstMate), "CGACCCATGT"),
        Read(ReadId("read1", MateNumber::kSecondMate), "GACCCATGTC"));

    regionAnalyzer.processMates(
        Read(ReadId("read2", MateNumber::kFirstMate), "CGACATGT"),
        Read(ReadId("read2", MateNumber::kSecondMate), "GACATGTC"));

    RegionFindings regionFindings = regionAnalyzer.analyze(sampleParams);

    std::unique_ptr<VariantFindings> repeatFindingsPtr(new RepeatFindings(
        CountTable({ { 1, 2 }, { 3, 2 } }), CountTable(), CountTable(), RepeatGenotype(1, { 1, 3 })));
    RegionFindings expectedFindings;
    expectedFindings.emplace(std::make_pair("repeat", std::move(repeatFindingsPtr)));

    ASSERT_EQ(expectedFindings, regionFindings);
}

TEST_P(AlignerTests, RegionAnalysis_ShortMultiUnitRepeat_Genotyped)
{
    //    const int32_t kmerLenForReadOrientation = 5;
    //    const int32_t kmerLenForAlignment = 4;
    //    const int32_t paddingLength = 1;
}

// INSTANTIATE_TEST_CASE_P(
//    AlignerTestsInst, AlignerTests, ::testing::Values(std::string("path-aligner"), std::string("dag-aligner")), );
