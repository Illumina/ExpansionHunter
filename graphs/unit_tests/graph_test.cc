//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
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

#include "graphs/graph.h"

#include <memory>

#include "gtest/gtest.h"

using std::set;

TEST(GraphConstruction, TypicalNodeCount_GraphConstructed) {
  Graph graph(3);
  ASSERT_EQ(3, graph.NumNodes());
}

TEST(NodeSequenceManipulation, TypicalSequence_SequenceSet) {
  Graph graph(3);
  graph.SetNodeSeq(1, "ATT");
  EXPECT_EQ("ATT", graph.NodeSeq(1));
}

TEST(NodeSequenceManipulation, NonexistingNode_ExceptionRaised) {
  Graph graph(3);
  EXPECT_ANY_THROW(graph.SetNodeSeq(4, "ATT"));
  EXPECT_ANY_THROW(graph.NodeSeq(4));
}

TEST(AddingEdgesToGraph, TypicalEdge_EdgeAdded) {
  Graph graph(3);
  graph.AddEdge(0, 1);
  graph.AddEdge(0, 0);
  EXPECT_TRUE(graph.HasEdge(0, 1));
  EXPECT_TRUE(graph.HasEdge(0, 0));
}

TEST(AddingEdgesToGraph, EdgeBreakingTopologicalOrder_ExceptionRaised) {
  Graph graph(3);
  EXPECT_ANY_THROW(graph.AddEdge(2, 1));
}

TEST(AddingEdgesToGraph, EdgesBetweenNonexistingNodes_ExceptionRaised) {
  Graph graph(4);
  EXPECT_ANY_THROW(graph.AddEdge(-1, 2));
  EXPECT_ANY_THROW(graph.AddEdge(1, 4));
  EXPECT_ANY_THROW(graph.AddEdge(4, 5));
}

TEST(CheckingEdgesInGraph, EdgesBetweenNonexistingNodes_ExceptionRaised) {
  Graph graph(4);
  EXPECT_ANY_THROW(graph.HasEdge(-1, 2));
  EXPECT_ANY_THROW(graph.HasEdge(1, 4));
  EXPECT_ANY_THROW(graph.HasEdge(4, 5));
}

TEST(GettingNodeNeighbors, TypicalNode_SuccessorsFound) {
  Graph graph(4);
  graph.AddEdge(0, 1);
  graph.AddEdge(0, 2);
  graph.AddEdge(0, 3);
  graph.AddEdge(2, 3);

  const set<int32_t> expected_successors = {1, 2, 3};
  ASSERT_EQ(expected_successors, graph.Successors(0));
  ASSERT_TRUE(graph.Successors(1).empty());
}

TEST(GettingNodeNeighbors, TypicalNode_PredecessorsFound) {
  Graph graph(4);
  graph.AddEdge(0, 1);
  graph.AddEdge(0, 2);
  graph.AddEdge(0, 3);
  graph.AddEdge(2, 3);

  const set<int32_t> expected_predecessors = {0, 2};
  ASSERT_EQ(expected_predecessors, graph.Predecessors(3));
}

TEST(GettingNodeNeighbors, NeighborsOfNonexistingNode_ExceptionRaised) {
  Graph graph(4);
  ASSERT_ANY_THROW(graph.Successors(4));
  ASSERT_ANY_THROW(graph.Predecessors(-1));
}

TEST(ResettingGraph, TypicalReset_NewGraphInitialized) {
  Graph graph(3);
  graph.AddEdge(0, 1);
  graph.Reset(2);
  EXPECT_EQ(2, graph.NumNodes());
  EXPECT_TRUE(graph.Successors(0).empty());
  EXPECT_TRUE(graph.Successors(1).empty());
}
