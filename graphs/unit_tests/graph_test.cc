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

struct AGraphWithFourNodes : ::testing::Test {
  void SetUp() {
    graph_ptr = std::unique_ptr<Graph>(new Graph(4));
    graph_ptr->AddEdge(0, 1);
    graph_ptr->AddEdge(0, 2);
    graph_ptr->AddEdge(0, 3);
    graph_ptr->AddEdge(2, 3);
  }
  std::unique_ptr<Graph> graph_ptr;
};

TEST(AGraph, IsConstructedByProvidingNodeCount) {
  Graph graph(3);
  ASSERT_EQ(3, graph.NumNodes());
}

TEST(AGraph, AllowsToAddEdgesAndCheckThatTheyExist) {
  Graph graph(3);
  graph.AddEdge(0, 1);
  ASSERT_TRUE(graph.HasEdge(0, 1));
}

TEST_F(AGraphWithFourNodes, ErrorsOutWhenAddingEdgesBreakingTopologicalOrder) {
  EXPECT_ANY_THROW(graph_ptr->AddEdge(2, 1));
}

TEST_F(AGraphWithFourNodes, ErrorsOutWhenAddingLoops) {
  EXPECT_ANY_THROW(graph_ptr->AddEdge(0, 0));
}

TEST_F(AGraphWithFourNodes, ErrorsOutWhenAddingEdgesWithNonexistingNodes) {
  EXPECT_ANY_THROW(graph_ptr->AddEdge(-1, 2));
  EXPECT_ANY_THROW(graph_ptr->AddEdge(1, 4));
  EXPECT_ANY_THROW(graph_ptr->AddEdge(4, 5));
}

TEST_F(AGraphWithFourNodes, ErrorsOutWhenCheckingForEdgesWithNonexistingNodes) {
  EXPECT_ANY_THROW(graph_ptr->HasEdge(-1, 2));
  EXPECT_ANY_THROW(graph_ptr->HasEdge(1, 4));
  EXPECT_ANY_THROW(graph_ptr->HasEdge(4, 5));
}

TEST_F(AGraphWithFourNodes, AllowsCheckingForSuccessorsOfNode) {
  const set<int32_t> expected_successors = {1, 2, 3};
  ASSERT_EQ(expected_successors, graph_ptr->Successors(0));
  ASSERT_TRUE(graph_ptr->Successors(1).empty());
}

TEST_F(AGraphWithFourNodes, AllowsCheckingForPredecessorsOfNode) {
  const set<int32_t> expected_predecessors = {0, 2};
  ASSERT_EQ(expected_predecessors, graph_ptr->Predecessors(3));
}

TEST_F(AGraphWithFourNodes, ErrorsOutWhenGettingNeighborsOfNonexistingNode) {
  Graph graph(4);
  ASSERT_ANY_THROW(graph_ptr->Successors(4));
  ASSERT_ANY_THROW(graph_ptr->Predecessors(-1));
}
