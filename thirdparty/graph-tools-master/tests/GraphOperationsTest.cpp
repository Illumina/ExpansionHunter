// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
