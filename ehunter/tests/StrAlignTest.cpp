//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
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

#include "genotyping/StrAlign.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using namespace ehunter;

TEST(GettingCompatibleAlignmentByClippingFromLeft, TypicalRead_CompatibleAlignmentFound)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    int strNodeId = 1;
    ConsistentAlignmentCalculator alignmentCalculator(strNodeId);

    GraphAlignment spanningAlign = decodeGraphAlignment(0, "0[5M2I1M]1[1M]1[1M]1[1M]2[1M1D2M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 2, 17, 0), alignmentCalculator.clipFromLeft(2, spanningAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kSpanning, 3, 36, 0), alignmentCalculator.clipFromLeft(3, spanningAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 22, 0), alignmentCalculator.clipFromLeft(4, spanningAlign));

    GraphAlignment rightFlankingAlign = decodeGraphAlignment(0, "1[1M]1[1M]1[1M]2[4M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 2, 30, 0), alignmentCalculator.clipFromLeft(2, rightFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 35, 0), alignmentCalculator.clipFromLeft(3, rightFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 35, 0), alignmentCalculator.clipFromLeft(4, rightFlankingAlign));

    GraphAlignment leftFlankingAlign = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]1[1M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 0, 0, 0), alignmentCalculator.clipFromLeft(0, leftFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 2, 10, 0), alignmentCalculator.clipFromLeft(2, leftFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 30, 0), alignmentCalculator.clipFromLeft(3, leftFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 30, 0), alignmentCalculator.clipFromLeft(4, leftFlankingAlign));

    GraphAlignment inRepeatAlign = decodeGraphAlignment(0, "1[1M]1[1M]1[1M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 0, 0, 0), alignmentCalculator.clipFromLeft(0, inRepeatAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 2, 10, 0), alignmentCalculator.clipFromLeft(2, inRepeatAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 3, 15, 0), alignmentCalculator.clipFromLeft(3, inRepeatAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 3, 15, 0), alignmentCalculator.clipFromLeft(4, inRepeatAlign));

    GraphAlignment alignInsideLeftFlank = decodeGraphAlignment(0, "0[6M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kOutside, 0, 30, 0), alignmentCalculator.clipFromLeft(2, alignInsideLeftFlank));

    GraphAlignment alignInsideRightFlank = decodeGraphAlignment(1, "2[5M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kOutside, 0, 25, 0), alignmentCalculator.clipFromLeft(2, alignInsideRightFlank));
}

TEST(GettingCompatibleAlignmentByClippingFromRight, TypicalRead_CompatibleAlignmentFound)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    ConsistentAlignmentCalculator alignmentCalculator(1);

    GraphAlignment spanningAlign = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]1[1M]2[4M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 2, 25, 0), alignmentCalculator.clipFromRight(2, spanningAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kSpanning, 3, 50, 0), alignmentCalculator.clipFromRight(3, spanningAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 30, 0), alignmentCalculator.clipFromRight(4, spanningAlign));

    GraphAlignment leftFlankingAlign = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]1[1M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 2, 25, 0), alignmentCalculator.clipFromRight(2, leftFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 30, 0), alignmentCalculator.clipFromRight(3, leftFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 30, 0), alignmentCalculator.clipFromRight(4, leftFlankingAlign));

    GraphAlignment rightFlankingAlign = decodeGraphAlignment(0, "1[1M]1[1M]1[1M]2[4M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 0, 0, 0), alignmentCalculator.clipFromRight(0, rightFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 2, 10, 0), alignmentCalculator.clipFromRight(2, rightFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 35, 0), alignmentCalculator.clipFromRight(3, rightFlankingAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 3, 35, 0), alignmentCalculator.clipFromRight(4, rightFlankingAlign));

    GraphAlignment inRepeatAlign = decodeGraphAlignment(0, "1[1M]1[1M]1[1M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 0, 0, 0), alignmentCalculator.clipFromRight(0, inRepeatAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 2, 10, 0), alignmentCalculator.clipFromRight(2, inRepeatAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 3, 15, 0), alignmentCalculator.clipFromRight(3, inRepeatAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kInRepeat, 3, 15, 0), alignmentCalculator.clipFromRight(4, inRepeatAlign));

    GraphAlignment alignInsideLeftFlank = decodeGraphAlignment(0, "0[6M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kOutside, 0, 30, 0), alignmentCalculator.clipFromRight(2, alignInsideLeftFlank));

    GraphAlignment alignInsideRightFlank = decodeGraphAlignment(1, "2[5M]", &graph);
    EXPECT_EQ(
        StrAlign(StrAlign::Type::kOutside, 0, 25, 0), alignmentCalculator.clipFromRight(2, alignInsideRightFlank));
}

TEST(GettingCompatibleAlignmentByRemovingStutter, TypicalRead_CompatibleAlignmentFound)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(AT)*ATGTCG"));
    ConsistentAlignmentCalculator alignmentCalculator(1);

    GraphAlignment flankingAlign = decodeGraphAlignment(0, "1[2M]1[2M]1[2M]2[4M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kOutside, 0, 0, 0), alignmentCalculator.removeStutter(3, flankingAlign));

    GraphAlignment spanningAlign = decodeGraphAlignment(3, "0[3M]1[2M]1[2M]1[2M]2[4M]", &graph);
    EXPECT_EQ(StrAlign(StrAlign::Type::kSpanning, 0, 0, 0), alignmentCalculator.removeStutter(0, spanningAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kSpanning, 1, 0, 0), alignmentCalculator.removeStutter(1, spanningAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kSpanning, 2, 19, 0), alignmentCalculator.removeStutter(2, spanningAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kSpanning, 3, 65, 0), alignmentCalculator.removeStutter(3, spanningAlign));
    EXPECT_EQ(StrAlign(StrAlign::Type::kSpanning, 4, 29, 0), alignmentCalculator.removeStutter(4, spanningAlign));
}

TEST(GettingCompatibleAlignment, TypicalRead_CompatibleAlignmentFound)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(AT)*ATGTCG"));
    ConsistentAlignmentCalculator alignmentCalculator(1);

    {
        GraphAlignment align = decodeGraphAlignment(3, "0[3M]1[2M]1[2M]1[2M]2[4M]", &graph);
        EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 2, 40, 0), alignmentCalculator.findConsistentAlignment(2, align));
    }

    {
        GraphAlignment align = decodeGraphAlignment(3, "0[3M]1[2M]1[2M]1[2M]2[2M]", &graph);
        EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 2, 35, 0), alignmentCalculator.findConsistentAlignment(2, align));
    }

    {
        GraphAlignment align = decodeGraphAlignment(0, "0[6M]1[2M]1[2M]1[2M]2[6M]", &graph);
        EXPECT_EQ(StrAlign(StrAlign::Type::kFlanking, 2, 50, 0), alignmentCalculator.findConsistentAlignment(2, align));
    }
}
