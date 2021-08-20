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

#include "io/RegionGraph.hh"

#include "gtest/gtest.h"

#include "graphcore/Graph.hh"

#include "io/GraphBlueprint.hh"
#include "locus/LocusSpecification.hh"

using graphtools::Graph;
using std::string;

using namespace ehunter;

TEST(ConstructingRepeatRegionGraphs, SingleUnitStr_GraphConstructed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    ASSERT_EQ(3u, graph.numNodes());
    ASSERT_EQ("ATTCGA", graph.nodeSeq(0));
    ASSERT_EQ("C", graph.nodeSeq(1));
    ASSERT_EQ("ATGTCG", graph.nodeSeq(2));

    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(1, 1));
    EXPECT_TRUE(graph.hasEdge(1, 2));
}

TEST(ConstructingRepeatRegionGraphs, MultiUnitStr_GraphConstructed)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AAAATT(AGG)*ATG(CG)*GGGGCC"));

    ASSERT_EQ(5u, graph.numNodes());
    EXPECT_EQ(8u, graph.numEdges());

    EXPECT_EQ("AAAATT", graph.nodeSeq(0));
    EXPECT_EQ("AGG", graph.nodeSeq(1));
    EXPECT_EQ("ATG", graph.nodeSeq(2));
    EXPECT_EQ("CG", graph.nodeSeq(3));
    EXPECT_EQ("GGGGCC", graph.nodeSeq(4));

    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(0, 2));

    EXPECT_TRUE(graph.hasEdge(1, 1));
    EXPECT_TRUE(graph.hasEdge(1, 2));

    EXPECT_TRUE(graph.hasEdge(2, 3));
    EXPECT_TRUE(graph.hasEdge(2, 4));

    EXPECT_TRUE(graph.hasEdge(3, 3));
    EXPECT_TRUE(graph.hasEdge(3, 4));
}
