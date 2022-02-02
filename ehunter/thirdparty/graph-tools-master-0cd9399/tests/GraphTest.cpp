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
