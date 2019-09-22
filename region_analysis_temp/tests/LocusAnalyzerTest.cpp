//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

#include "region_analysis/LocusAnalyzer.hh"

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
    initializeWorkflowContext(HeuristicParameters(1000, 20, true, GetParam(), 4, 1, 5));

    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    vector<GenomicRegion> referenceRegions = { GenomicRegion(1, 1, 2) };

    NodeToRegionAssociation dummyAssociation;
    GenotyperParameters params;
    LocusSpecification locusSpec("workflow", ChromType::kAutosome, referenceRegions, graph, dummyAssociation, params);
    VariantClassification classification(VariantType::kRepeat, VariantSubtype::kCommonRepeat);
    locusSpec.addVariantSpecification("repeat", classification, GenomicRegion(1, 1, 2), { 1 }, 1);

    graphtools::BlankAlignmentWriter blankAlignmentWriter;
    LocusAnalyzer locusAnalyzer(locusSpec, blankAlignmentWriter);

    locusAnalyzer.processMates(
        Read(ReadId("read1", MateNumber::kFirstMate), "CGACCCATGT", true),
        Read(ReadId("read1", MateNumber::kSecondMate), "GACCCATGTC", true), RegionType::kTarget);

    locusAnalyzer.processMates(
        Read(ReadId("read2", MateNumber::kFirstMate), "CGACATGT", true),
        Read(ReadId("read2", MateNumber::kSecondMate), "GACATGTC", true), RegionType::kTarget);

    LocusFindings locusFindings = locusAnalyzer.analyze(Sex::kFemale);

    std::unique_ptr<VariantFindings> repeatFindingsPtr(new RepeatFindings(
        CountTable({ { 1, 2 }, { 3, 2 } }), CountTable(), CountTable(), RepeatGenotype(1, { 1, 3 })));
    LocusFindings expectedFindings;
    expectedFindings.findingsForEachVariant.emplace("repeat", std::move(repeatFindingsPtr));

    ASSERT_EQ(expectedFindings.findingsForEachVariant, locusFindings.findingsForEachVariant);
}

TEST_P(AlignerTests, RegionAnalysis_ShortMultiUnitRepeat_Genotyped)
{
    //    const int32_t kmerLenForReadOrientation = 5;
    //    const int32_t kmerLenForAlignment = 4;
    //    const int32_t paddingLength = 1;
}

// INSTANTIATE_TEST_CASE_P(
//    AlignerTestsInst, AlignerTests, ::testing::Values(std::string("path-aligner"), std::string("dag-aligner")), );
