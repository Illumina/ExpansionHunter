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

#include "gtest/gtest.h"

TEST(AGraph, IsConstructedByProvidingNodeCount) {
  Graph graph(3);
  ASSERT_EQ(3, graph.NumNodes());
}

TEST(AGraph, AllowsToAddEdgesAndCheckThatTheyExist) {
  Graph graph(3);
  graph.AddEdge(0, 1);
  ASSERT_TRUE(graph.HasEdge(0, 1));
}

TEST(AGraph, ErrorsOutWhenAddingEdgesBreakingTopologicalOrderAndLoops) {
  Graph graph(3);
  EXPECT_ANY_THROW(graph.AddEdge(0, 0));
  EXPECT_ANY_THROW(graph.AddEdge(2, 1));
}

TEST(AGraph, ErrorsOutWhenAddingEdgesWithNonexistingNodes) {
  Graph graph(4);
  EXPECT_ANY_THROW(graph.AddEdge(-1, 2));
  EXPECT_ANY_THROW(graph.AddEdge(1, 4));
  EXPECT_ANY_THROW(graph.AddEdge(4, 5));
}

TEST(AGraph, ErrorsOutWhenCheckingForEdgesWithNonexistingNodes) {
  Graph graph(4);
  EXPECT_ANY_THROW(graph.HasEdge(-1, 2));
  EXPECT_ANY_THROW(graph.HasEdge(1, 4));
  EXPECT_ANY_THROW(graph.HasEdge(4, 5));
}