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
  const string del = "CCCC";
  const string right_flank = "GGGCC";
  const Graph graph = makeDeletionGraph(left_flank, del, right_flank);

  EXPECT_EQ(3, graph.NumNodes());
  EXPECT_EQ(left_flank, graph.NodeSeq(0));
  EXPECT_EQ(del, graph.NodeSeq(1));
  EXPECT_EQ(right_flank, graph.NodeSeq(2));
  EXPECT_TRUE(graph.HasEdge(0, 1));
  EXPECT_TRUE(graph.HasEdge(0, 2));
  EXPECT_TRUE(graph.HasEdge(1, 2));
}