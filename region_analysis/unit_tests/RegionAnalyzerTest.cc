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

#include "region_spec/RegionBlueprint.hh"
#include "region_spec/RegionSpec.hh"

using graphtools::Graph;
using reads::Read;
using std::string;
using std::vector;

class AlignerTests : public ::testing::TestWithParam<std::string> {};

// TODO: Throw an error if there are no valid extensions?
TEST_P(AlignerTests, RegionAnalysis_ShortSingleUnitRepeat_Genotyped)
//TEST(DISABLED_RegionAnalysis, ShortSingleUnitRepeat_Genotyped)
{
    const vector<string> repeatIds;

    RegionBlueprint blueprint(
        "ATTCGA", "C", "ATGTCG", { "repeat" }, { Region("chr1:1-2") }, { RegionBlueprintComponent::Rarity::kCommon });

    RegionSpec regionSpec("region", blueprint, AlleleCount::kTwo, Region("chr1:1-2"));

    const int readLength = 10;
    const double haplotypeDepth = 5.0;

    GraphAlignmentHeuristicsParameters alignmentParams(4, 1, 5);

    RegionAnalyzer regionAnalyzer(regionSpec, haplotypeDepth, readLength, std::cerr, GetParam(), alignmentParams);

    regionAnalyzer.processMates(Read("read1/1", "CGACCCATGT"), Read("read1/2", "GACCCATGTC"));
    regionAnalyzer.processMates(Read("read2/1", "CGACATGT"), Read("read2/2", "GACATGTC"));

    RegionFindings regionFindings = regionAnalyzer.genotype();

    RepeatFindings repeatFindings(CountTable(), CountTable({ { 1, 2 }, { 3, 2 } }), RepeatGenotype(1, { 1, 3 }));
    RegionFindings expectedFindings = { { "repeat", repeatFindings } };

    ASSERT_EQ(expectedFindings, regionFindings);
}

TEST_P(AlignerTests, RegionAnalysis_ShortMultiUnitRepeat_Genotyped)
{
//    const int32_t kmerLenForReadOrientation = 5;
//    const int32_t kmerLenForAlignment = 4;
//    const int32_t paddingLength = 1;
}

//INSTANTIATE_TEST_CASE_P(
//    AlignerTestsInst, AlignerTests, ::testing::Values(std::string("path-aligner"), std::string("dag-aligner")), );
