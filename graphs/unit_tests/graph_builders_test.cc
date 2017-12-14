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

#include "graphs/graph_builders.h"

#include "gtest/gtest.h"

using std::string;

TEST(DeletionGraph, IsCreatedFromNodeSequences) {
  const string left_flank = "AATT";
  const string deletion = "CCCC";
  const string right_flank = "GGGCC";
  const Graph graph = makeDeletionGraph(left_flank, deletion, right_flank);

  EXPECT_EQ(3, graph.NumNodes());
  EXPECT_EQ(left_flank, graph.NodeSeq(0));
  EXPECT_EQ(deletion, graph.NodeSeq(1));
  EXPECT_EQ(right_flank, graph.NodeSeq(2));
  EXPECT_TRUE(graph.HasEdge(0, 1));
  EXPECT_TRUE(graph.HasEdge(0, 2));
  EXPECT_TRUE(graph.HasEdge(1, 2));
}

TEST(SwapGraph, IsCreatedFromNodeSequences) {
  const string left_flank = "AATT";
  const string deletion = "CCCC";
  const string insertion = "TTTT";
  const string right_flank = "GGGCC";
  const Graph graph =
      makeSwapGraph(left_flank, deletion, insertion, right_flank);

  EXPECT_EQ(4, graph.NumNodes());
  EXPECT_EQ(left_flank, graph.NodeSeq(0));
  EXPECT_EQ(deletion, graph.NodeSeq(1));
  EXPECT_EQ(insertion, graph.NodeSeq(2));
  EXPECT_EQ(right_flank, graph.NodeSeq(3));
  EXPECT_TRUE(graph.HasEdge(0, 1));
  EXPECT_TRUE(graph.HasEdge(0, 2));
  EXPECT_TRUE(graph.HasEdge(1, 3));
  EXPECT_TRUE(graph.HasEdge(2, 3));
}

TEST(DoubleSwapGraph, IsCreatedFromNodeSequences) {
  const string left_flank = "AATT";
  const string deletion1 = "CCCC";
  const string insertion1 = "TTTT";
  const string middle = "CCCC";
  const string deletion2 = "AAAA";
  const string insertion2 = "GGGG";
  const string right_flank = "GGGCC";
  const Graph graph =
      makeDoubleSwapGraph(left_flank, deletion1, insertion1, middle, deletion2,
                          insertion2, right_flank);

  EXPECT_EQ(7, graph.NumNodes());
  EXPECT_EQ(left_flank, graph.NodeSeq(0));
  EXPECT_EQ(deletion1, graph.NodeSeq(1));
  EXPECT_EQ(insertion1, graph.NodeSeq(2));
  EXPECT_EQ(middle, graph.NodeSeq(3));
  EXPECT_EQ(deletion2, graph.NodeSeq(4));
  EXPECT_EQ(insertion2, graph.NodeSeq(5));
  EXPECT_EQ(right_flank, graph.NodeSeq(6));
  EXPECT_TRUE(graph.HasEdge(0, 1));
  EXPECT_TRUE(graph.HasEdge(0, 2));
  EXPECT_TRUE(graph.HasEdge(1, 3));
  EXPECT_TRUE(graph.HasEdge(2, 3));
  EXPECT_TRUE(graph.HasEdge(3, 4));
  EXPECT_TRUE(graph.HasEdge(3, 5));
  EXPECT_TRUE(graph.HasEdge(4, 6));
  EXPECT_TRUE(graph.HasEdge(5, 6));
}

TEST(ConstructionOfLooplessStrGraph, TypicalParameters_GraphConstructed) {
  const string left_flank = "AATT";
  const string repeat_unit = "CGG";
  const string right_flank = "ATTT";
  const int32_t read_len = 10;
  const Graph graph =
      makeLooplessStrGraph(read_len, left_flank, repeat_unit, right_flank);

  ASSERT_EQ(6, graph.NumNodes());
  EXPECT_EQ(left_flank, graph.NodeSeq(0));
  EXPECT_EQ(repeat_unit, graph.NodeSeq(1));
  EXPECT_EQ(repeat_unit, graph.NodeSeq(2));
  EXPECT_EQ(repeat_unit, graph.NodeSeq(3));
  EXPECT_EQ(repeat_unit, graph.NodeSeq(4));
  EXPECT_EQ(right_flank, graph.NodeSeq(5));

  EXPECT_TRUE(graph.HasEdge(0, 5));

  EXPECT_TRUE(graph.HasEdge(0, 1));
  EXPECT_TRUE(graph.HasEdge(1, 5));

  EXPECT_TRUE(graph.HasEdge(1, 2));
  EXPECT_TRUE(graph.HasEdge(2, 5));

  EXPECT_TRUE(graph.HasEdge(2, 3));
  EXPECT_TRUE(graph.HasEdge(3, 5));

  EXPECT_TRUE(graph.HasEdge(3, 4));
  EXPECT_TRUE(graph.HasEdge(4, 5));
}

TEST(ConstructionOfStrGraph, TypicalParameters_GraphConstructed) {
  const string left_flank = "AATT";
  const string repeat_unit = "CGG";
  const string right_flank = "ATTT";
  const Graph graph = makeStrGraph(left_flank, repeat_unit, right_flank);

  ASSERT_EQ(3, graph.NumNodes());
  EXPECT_EQ(left_flank, graph.NodeSeq(0));
  EXPECT_EQ(repeat_unit, graph.NodeSeq(1));
  EXPECT_EQ(right_flank, graph.NodeSeq(2));

  EXPECT_TRUE(graph.HasEdge(0, 1));
  EXPECT_TRUE(graph.HasEdge(0, 2));
  EXPECT_TRUE(graph.HasEdge(1, 1));
  EXPECT_TRUE(graph.HasEdge(1, 2));
}
