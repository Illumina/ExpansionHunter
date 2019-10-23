//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "locus_spec/LocusSpecification.hh"

#include "gtest/gtest.h"

using std::string;
using std::vector;

using namespace ehunter;

/*
TEST(CreatingRegionSpecs, SingeUnitRepeatSpec_Parsed)
{
    const string locusId = "region1";
    const vector<string> repeatIds;
    RegionSpec locusSpec(locusId, repeatIds, Region("chr1:1-10"), "GCC", "AT", "GCC", "TA");

    vector<RegionComponent> expectedComponents = { RegionComponent("LF", "AT", RegionComponentType::kFlank),
                                                   RegionComponent("Repeat0", "GCC", RegionComponentType::kRepeat),
                                                   RegionComponent("RF", "TA", RegionComponentType::kFlank) };

    ASSERT_EQ(expectedComponents, locusSpec.components());
}

TEST(CreatingRegionSpecs, MultiunitFromatToSpecifySingleUnitRepeat_Parsed)
{
    const string locusId = "region1";
    const vector<string> repeatIds;
    RegionSpec locusSpec(locusId, repeatIds, Region("chr1:1-10"), "(AGG)", "AT", "AGG", "TA");

    vector<RegionComponent> expectedComponents = { RegionComponent("LF", "AT", RegionComponentType::kFlank),
                                                   RegionComponent("Repeat0", "AGG", RegionComponentType::kRepeat),
                                                   RegionComponent("RF", "TA", RegionComponentType::kFlank) };

    ASSERT_EQ(expectedComponents, locusSpec.components());
}

TEST(CreatingRegionSpecs, TypicalMultiunitRepeatSpecs_Parsed)
{
    const string locusId = "region1";
    const vector<string> repeatIds;

    {
        RegionSpec locusSpec(locusId, repeatIds, Region("chr1:1-10"), "(AGG)(CG)", "AT", "AGGCG", "TA");

        vector<RegionComponent> expectedComponents = { RegionComponent("LF", "AT", RegionComponentType::kFlank),
                                                       RegionComponent("Repeat0", "AGG", RegionComponentType::kRepeat),
                                                       RegionComponent("Repeat1", "CG", RegionComponentType::kRepeat),
                                                       RegionComponent("RF", "TA", RegionComponentType::kFlank) };

        EXPECT_EQ(expectedComponents, locusSpec.components());
    }

    {
        RegionSpec locusSpec(locusId, repeatIds, Region("chr1:1-10"), "(AGG)ATG(CG)", "CC", "AGGATGCG", "GG");

        vector<RegionComponent> expectedComponents = { RegionComponent("LF", "CC", RegionComponentType::kFlank),
                                                       RegionComponent("Repeat0", "AGG", RegionComponentType::kRepeat),
                                                       RegionComponent("", "ATG", RegionComponentType::kInterruption),
                                                       RegionComponent("Repeat1", "CG", RegionComponentType::kRepeat),
                                                       RegionComponent("RF", "GG", RegionComponentType::kFlank) };

        EXPECT_EQ(expectedComponents, locusSpec.components());
    }
}
*/
