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

#include "genotyping/AlignMatrix.hh"

#include "gmock/gmock.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using std::vector;

using namespace ehunter;
using namespace strgt;

TEST(CreatingAlignmentMatrix, LikelihoodMatrixInitialization_Initialized)
{
    AlignMatrix alignMatrix(0);
    ASSERT_EQ(0, alignMatrix.numReads());
}

TEST(ReadLikelihoods, SpanningRead_Computed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    int strNodeId = 1;
    AlignMatrix alignMatrix(strNodeId);
    GraphAlignment read = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]2[4M]", &graph);
    GraphAlignment mate = decodeGraphAlignment(0, "0[6M]", &graph);

    alignMatrix.add(read, mate);
    ASSERT_EQ(alignMatrix.getAlign(0, 0), StrAlign('F', 0, 20, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 1), StrAlign('F', 1, 25, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 2), StrAlign('S', 2, 45, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 3), StrAlign('F', 2, 30, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 4), StrAlign('F', 2, 30, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 5), StrAlign('F', 2, 30, 0));
}

TEST(ReadLikelihoods, LeftFlankingRead_Computed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    AlignMatrix alignMatrix(1);
    GraphAlignment read = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]", &graph);
    GraphAlignment mate = decodeGraphAlignment(0, "0[6M]", &graph);

    alignMatrix.add(read, mate);
    ASSERT_EQ(alignMatrix.getAlign(0, 0), StrAlign('F', 0, 15, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 1), StrAlign('F', 1, 20, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 2), StrAlign('F', 2, 25, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 3), StrAlign('F', 2, 25, 0));
}

TEST(ReadLikelihoods, RightFlankingRead_Computed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    AlignMatrix alignMatrix(1);
    GraphAlignment read = decodeGraphAlignment(0, "1[1M]1[1M]1[1M]2[4M]", &graph);
    GraphAlignment mate = decodeGraphAlignment(0, "0[6M]", &graph);

    alignMatrix.add(read, mate);
    ASSERT_EQ(alignMatrix.getAlign(0, 0), StrAlign('F', 0, 20, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 1), StrAlign('F', 1, 25, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 2), StrAlign('F', 2, 30, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 3), StrAlign('F', 3, 35, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 4), StrAlign('F', 3, 35, 0));
}

TEST(ReadLikelihoods, InRepeatRead_Computed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    AlignMatrix alignMatrix(1);
    GraphAlignment read = decodeGraphAlignment(0, "1[1M]1[1M]1[1M]1[1M]", &graph);
    GraphAlignment mate = decodeGraphAlignment(0, "0[6M]", &graph);

    alignMatrix.add(read, mate);
    ASSERT_EQ(alignMatrix.getAlign(0, 0), StrAlign('I', 0, 0, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 1), StrAlign('I', 1, 5, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 2), StrAlign('I', 2, 10, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 3), StrAlign('I', 3, 15, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 4), StrAlign('I', 4, 20, 0));
    ASSERT_EQ(alignMatrix.getAlign(0, 5), StrAlign('I', 4, 20, 0));
}

TEST(AddingIrrPairs, NoOtherIrrsPresent_IrrPairsIgnored)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    int strNodeId = 1;
    AlignMatrix alignMatrix(strNodeId);
    GraphAlignment mate = decodeGraphAlignment(0, "0[5M]", &graph);
    alignMatrix.add(decodeGraphAlignment(0, "1[1M]1[1M]1[1M]", &graph), mate);

    int maxMotifsInRead = 5;
    int numIrrPairs = 5;
    addIrrPairsIfPossibleExpansion(maxMotifsInRead, alignMatrix, numIrrPairs);
    ASSERT_EQ(alignMatrix.numReads(), 2);
}

TEST(AddingIrrPairs, OtherIrrsPresent_IrrPairsAdded)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    int strNodeId = 1;
    AlignMatrix alignMatrix(strNodeId);

    GraphAlignment mate = decodeGraphAlignment(0, "0[5M]", &graph);
    alignMatrix.add(decodeGraphAlignment(0, "1[1M]1[1M]1[1M]1[1M]1[1M]", &graph), mate);

    int maxMotifsInRead = 6;
    int numIrrPairs = 2;
    addIrrPairsIfPossibleExpansion(maxMotifsInRead, alignMatrix, numIrrPairs);
    ASSERT_EQ(alignMatrix.numReads(), 6);
}
