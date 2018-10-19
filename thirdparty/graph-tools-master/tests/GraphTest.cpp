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

#include "graphcore/Graph.hh"

#include "gtest/gtest.h"

using std::set;
using std::string;
using std::vector;

using namespace graphtools;

TEST(GraphConstruction, TypicalNodeCount_GraphConstructed)
{
    Graph graph(3);
    ASSERT_EQ(3ul, graph.numNodes());
}

TEST(NodeNameManipulation, TypicalNode_NodeNameSet)
{
    Graph graph(3);
    graph.setNodeName(1, "LF");
    ASSERT_EQ("LF", graph.nodeName(1));
}

TEST(NodeNameManipulation, NonexistingNode_ExceptionRaised)
{
    Graph graph(1);
    ASSERT_ANY_THROW(graph.setNodeName(1, "LF"));
    ASSERT_ANY_THROW(graph.nodeName(1));
}

TEST(NodeSequenceManipulation, TypicalSequence_SequenceSet)
{
    Graph graph(3);
    graph.setNodeSeq(1, "ATT");
    EXPECT_EQ("ATT", graph.nodeSeq(1));
}

TEST(NodeSequenceManipulation, DegenerateSequence_SequenceExpansionObtained)
{
    Graph graph(3);
    graph.setNodeSeq(1, "WC");
    vector<string> expected_expansion = { "AC", "TC" };
    EXPECT_EQ(expected_expansion, graph.nodeSeqExpansion(1));
}

TEST(NodeSequenceManipulation, NonexistingNode_ExceptionRaised)
{
    Graph graph(3);
    EXPECT_ANY_THROW(graph.setNodeSeq(4, "ATT"));
    EXPECT_ANY_THROW(graph.nodeSeq(4));
    EXPECT_ANY_THROW(graph.nodeSeqExpansion(4));
}

TEST(NodeSequenceManipulation, EmptySequence_ExceptionRaised)
{
    Graph graph(3);
    EXPECT_ANY_THROW(graph.setNodeSeq(1, ""));
}

TEST(AddingEdges, TypicalEdge_EdgeAdded)
{
    Graph graph(3);
    graph.addEdge(0, 1);
    graph.addEdge(0, 0);
    EXPECT_TRUE(graph.hasEdge(0, 1));
    EXPECT_TRUE(graph.hasEdge(0, 0));
}

TEST(AddingEdges, EdgeBreakingTopologicalOrder_ExceptionRaised)
{
    Graph graph(3);
    EXPECT_ANY_THROW(graph.addEdge(2, 1));
}

TEST(AddingEdges, EdgesBetweenNonexistingNodes_ExceptionRaised)
{
    Graph graph(4);
    EXPECT_ANY_THROW(graph.addEdge(-1, 2));
    EXPECT_ANY_THROW(graph.addEdge(1, 4));
    EXPECT_ANY_THROW(graph.addEdge(4, 5));
}

TEST(AddingEdges, EdgesThatAlreadyExist_ExceptionRaised)
{
    Graph graph(4);
    graph.addEdge(1, 2);
    EXPECT_ANY_THROW(graph.addEdge(1, 2));
}

TEST(CheckingIfEdgesExist, EdgesBetweenNonexistingNodes_ExceptionRaised)
{
    Graph graph(4);
    EXPECT_ANY_THROW(graph.hasEdge(-1, 2));
    EXPECT_ANY_THROW(graph.hasEdge(1, 4));
    EXPECT_ANY_THROW(graph.hasEdge(4, 5));
}

TEST(EdgeLabelManipulation, TypicalEdges_EdgesLabeled)
{
    Graph graph(4);
    graph.addEdge(0, 2);
    graph.addLabelToEdge(0, 2, "ref");
    graph.addLabelToEdge(0, 2, "alt");
    Labels expected_labels = { "ref", "alt" };
    ASSERT_EQ(expected_labels, graph.edgeLabels(0, 2));
}

TEST(EdgeLabelManipulation, NonexistingEdges_ExceptionRaised)
{
    Graph graph(4);
    EXPECT_ANY_THROW(graph.addLabelToEdge(0, 1, "ref"));
    EXPECT_ANY_THROW(graph.addLabelToEdge(0, 4, "ref"));
    EXPECT_ANY_THROW(graph.edgeLabels(0, 1));
}

TEST(GettingNodeNeighbors, TypicalNode_SuccessorsFound)
{
    Graph graph(4);
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(0, 3);
    graph.addEdge(2, 3);

    const set<NodeId> expected_successors = { 1, 2, 3 };
    ASSERT_EQ(expected_successors, graph.successors(0));
    ASSERT_TRUE(graph.successors(1).empty());
}

TEST(GettingNodeNeighbors, LoopAtNode_SuccessorsFound)
{
    Graph graph(4);
    graph.addEdge(0, 0);
    graph.addEdge(0, 1);

    const set<NodeId> expected_predecessors = { 0, 1 };
    ASSERT_EQ(expected_predecessors, graph.successors(0));
}

TEST(GettingNodeNeighbors, TypicalNode_PredecessorsFound)
{
    Graph graph(4);
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(0, 3);
    graph.addEdge(2, 3);

    const set<NodeId> expected_predecessors = { 0, 2 };
    ASSERT_EQ(expected_predecessors, graph.predecessors(3));
}

TEST(GettingNodeNeighbors, NeighborsOfNonexistingNode_ExceptionRaised)
{
    Graph graph(4);
    EXPECT_ANY_THROW(graph.successors(4));
    EXPECT_ANY_THROW(graph.predecessors(-1));
}
