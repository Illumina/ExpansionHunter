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

#include "input/RegionGraph.hh"

#include "gtest/gtest.h"

#include "graphcore/Graph.hh"

#include "input/GraphBlueprint.hh"
#include "region_spec/LocusSpecification.hh"

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
