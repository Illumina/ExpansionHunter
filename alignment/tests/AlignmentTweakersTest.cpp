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

#include "alignment/AlignmentTweakers.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"

#include "input/GraphBlueprint.hh"
#include "input/RegionGraph.hh"
#include "locus_spec/LocusSpecification.hh"

using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using std::string;

using namespace ehunter;

TEST(ShrinkingAlignmentPrefix, AlignmentWithUncertainPrefix_Shrank)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("CATGGTGA(A)*(GAA)*TAACTACT"));

    //                    --22222233333
    const string query = "TTGAAGAATAACT";

    {
        GraphAlignment alignment = decodeGraphAlignment(0, "2[2S3M]2[3M]3[5M]", &graph);
        shrinkUncertainPrefix(4, query, alignment);

        GraphAlignment expectedAlignment = decodeGraphAlignment(0, "2[5S3M]3[5M]", &graph);

        EXPECT_EQ(expectedAlignment, alignment);
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(0, "2[2S3M]2[3M]3[5M]", &graph);
        shrinkUncertainPrefix(8, query, alignment);

        GraphAlignment expectedAlignment = decodeGraphAlignment(0, "3[8S5M]", &graph);

        EXPECT_EQ(expectedAlignment, alignment);
    }
}

TEST(DISABLED_ShrinkingAlignmentSuffix, AlignmentWithUncertainSuffix_Shrank)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("CATGGTGA(A)*(GAA)*TAACTACT"));

    //                    0000011333--
    const string query = "GGTGAAATAAGG";
    GraphAlignment alignment = decodeGraphAlignment(3, "0[5M]1[1M]1[1M]3[3M2S]", &graph);

    shrinkUncertainSuffix(4, query, alignment);

    const GraphAlignment expectedAlignment = decodeGraphAlignment(3, "0[5M]1[1M]1[1M5S]", &graph);

    EXPECT_EQ(expectedAlignment, alignment);
}
