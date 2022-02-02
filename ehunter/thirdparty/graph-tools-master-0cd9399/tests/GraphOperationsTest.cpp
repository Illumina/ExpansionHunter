//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
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

#include "graphcore/GraphOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"

#include "gtest/gtest.h"

using std::set;
using std::string;
using std::vector;

using namespace graphtools;

TEST(GraphReversal, SwapGraph_Reversed)
{
    Graph graph = makeSwapGraph("CCCC", "AAAA", "GGGG", "TTTT");
    ASSERT_EQ(4ul, graph.numNodes());
    auto reversed_graph = reverseGraph(graph, false);
    ASSERT_EQ(4ul, reversed_graph.numNodes());
    ASSERT_EQ("TTTT", reversed_graph.nodeSeq(0));
    ASSERT_EQ("GGGG", reversed_graph.nodeSeq(1));
    ASSERT_EQ("AAAA", reversed_graph.nodeSeq(2));
    ASSERT_EQ("CCCC", reversed_graph.nodeSeq(3));
    ASSERT_EQ((set<NodeId>{ 1, 2 }), reversed_graph.successors(0));
    ASSERT_EQ((set<NodeId>{ 3 }), reversed_graph.successors(1));
    ASSERT_EQ((set<NodeId>{ 3 }), reversed_graph.successors(2));
    ASSERT_TRUE(reversed_graph.successors(3).empty());
}

TEST(GraphReversal, SwapGraph_SequenceReversed)
{
    Graph graph = makeSwapGraph("ACCC", "ATAA", "GGTG", "TTTA");
    ASSERT_EQ(4ul, graph.numNodes());
    auto reversed_graph = reverseGraph(graph, false);
    ASSERT_EQ(4ul, reversed_graph.numNodes());
    ASSERT_EQ("ATTT", reversed_graph.nodeSeq(0));
    ASSERT_EQ("GTGG", reversed_graph.nodeSeq(1));
    ASSERT_EQ("AATA", reversed_graph.nodeSeq(2));
    ASSERT_EQ("CCCA", reversed_graph.nodeSeq(3));
    ASSERT_EQ((set<NodeId>{ 1, 2 }), reversed_graph.successors(0));
    ASSERT_EQ((set<NodeId>{ 3 }), reversed_graph.successors(1));
    ASSERT_EQ((set<NodeId>{ 3 }), reversed_graph.successors(2));
    ASSERT_TRUE(reversed_graph.successors(3).empty());
}

TEST(GraphReversal, SwapGraph_SequenceReverseComplemented)
{
    Graph graph = makeSwapGraph("ACCC", "ATAA", "GGTG", "TTTA");
    ASSERT_EQ(4ul, graph.numNodes());
    auto reversed_graph = reverseGraph(graph, true);
    ASSERT_EQ(4ul, reversed_graph.numNodes());
    ASSERT_EQ("TAAA", reversed_graph.nodeSeq(0));
    ASSERT_EQ("CACC", reversed_graph.nodeSeq(1));
    ASSERT_EQ("TTAT", reversed_graph.nodeSeq(2));
    ASSERT_EQ("GGGT", reversed_graph.nodeSeq(3));
    ASSERT_EQ((set<NodeId>{ 1, 2 }), reversed_graph.successors(0));
    ASSERT_EQ((set<NodeId>{ 3 }), reversed_graph.successors(1));
    ASSERT_EQ((set<NodeId>{ 3 }), reversed_graph.successors(2));
    ASSERT_TRUE(reversed_graph.successors(3).empty());
}
