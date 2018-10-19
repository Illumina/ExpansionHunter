//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "region_spec/RegionSpec.hh"

#include "gtest/gtest.h"

using std::string;
using std::vector;

/*
TEST(CreatingRegionSpecs, SingeUnitRepeatSpec_Parsed)
{
    const string regionId = "region1";
    const vector<string> repeatIds;
    RegionSpec regionSpec(regionId, repeatIds, Region("chr1:1-10"), "GCC", "AT", "GCC", "TA");

    vector<RegionComponent> expectedComponents = { RegionComponent("LF", "AT", RegionComponentType::kFlank),
                                                   RegionComponent("Repeat0", "GCC", RegionComponentType::kRepeat),
                                                   RegionComponent("RF", "TA", RegionComponentType::kFlank) };

    ASSERT_EQ(expectedComponents, regionSpec.components());
}

TEST(CreatingRegionSpecs, MultiunitFromatToSpecifySingleUnitRepeat_Parsed)
{
    const string regionId = "region1";
    const vector<string> repeatIds;
    RegionSpec regionSpec(regionId, repeatIds, Region("chr1:1-10"), "(AGG)", "AT", "AGG", "TA");

    vector<RegionComponent> expectedComponents = { RegionComponent("LF", "AT", RegionComponentType::kFlank),
                                                   RegionComponent("Repeat0", "AGG", RegionComponentType::kRepeat),
                                                   RegionComponent("RF", "TA", RegionComponentType::kFlank) };

    ASSERT_EQ(expectedComponents, regionSpec.components());
}

TEST(CreatingRegionSpecs, TypicalMultiunitRepeatSpecs_Parsed)
{
    const string regionId = "region1";
    const vector<string> repeatIds;

    {
        RegionSpec regionSpec(regionId, repeatIds, Region("chr1:1-10"), "(AGG)(CG)", "AT", "AGGCG", "TA");

        vector<RegionComponent> expectedComponents = { RegionComponent("LF", "AT", RegionComponentType::kFlank),
                                                       RegionComponent("Repeat0", "AGG", RegionComponentType::kRepeat),
                                                       RegionComponent("Repeat1", "CG", RegionComponentType::kRepeat),
                                                       RegionComponent("RF", "TA", RegionComponentType::kFlank) };

        EXPECT_EQ(expectedComponents, regionSpec.components());
    }

    {
        RegionSpec regionSpec(regionId, repeatIds, Region("chr1:1-10"), "(AGG)ATG(CG)", "CC", "AGGATGCG", "GG");

        vector<RegionComponent> expectedComponents = { RegionComponent("LF", "CC", RegionComponentType::kFlank),
                                                       RegionComponent("Repeat0", "AGG", RegionComponentType::kRepeat),
                                                       RegionComponent("", "ATG", RegionComponentType::kInterruption),
                                                       RegionComponent("Repeat1", "CG", RegionComponentType::kRepeat),
                                                       RegionComponent("RF", "GG", RegionComponentType::kFlank) };

        EXPECT_EQ(expectedComponents, regionSpec.components());
    }
}
*/