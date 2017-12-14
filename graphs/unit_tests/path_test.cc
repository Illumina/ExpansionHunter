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

#include "graphs/path.h"

#include "gtest/gtest.h"

#include "graphs/graph.h"
#include "graphs/graph_builders.h"

using std::list;
using std::string;
using std::vector;

TEST(GettingPathSequence, TypicalPathOnDeletionGraph_SequenceReturned) {
  Graph graph = makeDeletionGraph("AAAACC", "TTTGG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  {
    GraphPath path(graph_ptr, 3, {0}, 3);
    EXPECT_EQ("A", path.seq());
  }

  {
    GraphPath path(graph_ptr, 3, {1}, 4);
    EXPECT_EQ("GG", path.seq());
  }

  {
    GraphPath path(graph_ptr, 3, {0, 1, 2}, 1);
    EXPECT_EQ("ACCTTTGGAT", path.seq());
  }
}

TEST(GettingPathSequence, TypicalPathOnStrGraph_SequenceReturned) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  GraphPath path(graph_ptr, 1, {0, 1, 1, 2}, 0);
  EXPECT_EQ("TTATATC", path.seq());
}

TEST(GettingLengthOfPathOverEachNode, TypicalPathOnStrGraph_LengthReturned) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 2, {0, 1, 1}, 0);

  EXPECT_EQ((size_t)1, path.GetOverlapWithNodeByIndex(0));
  EXPECT_EQ((size_t)2, path.GetOverlapWithNodeByIndex(1));
  EXPECT_EQ((size_t)1, path.GetOverlapWithNodeByIndex(2));
}

/*
TEST(GettingPathLength, TypicalPathOnStrGraph_LengthReturned) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  {
    GraphPath path(graph_ptr, 2, {0}, 2);
    EXPECT_EQ((size_t)1, path.length());
  }

  {
    GraphPath path(graph_ptr, 0, {1}, 1);
    EXPECT_EQ((size_t)2, path.length());
  }

  {
    GraphPath path(graph_ptr, 2, {0, 1, 1}, 0);
    EXPECT_EQ((size_t)4, path.length());
  }
}*/

/*
TEST_F(DeletionGraph, PathReturnsSequenceOfNodeItOverlaps) {
  GraphPath path(graph_ptr, 3, {0, 1, 2}, 1);
  EXPECT_EQ("ACC", path.seqOnNode(0));
  EXPECT_EQ("TTTGG", path.seqOnNode(1));
  EXPECT_EQ("AT", path.seqOnNode(2));
}

TEST_F(DeletionGraph, WellFormedPathIsValid) {
  GraphPath path(graph_ptr, 3, {0, 1, 2}, 1);
  ASSERT_TRUE(path.isValid());
}

TEST_F(DeletionGraph, PathStartingOutsideOfNodeSequenceIsInvalid) {
  GraphPath path(graph_ptr, 6, {0, 1, 2}, 1);
  ASSERT_FALSE(path.isValid());
}

TEST_F(DeletionGraph, PathEndingOutsideOfNodeSequenceIsInvalid) {
  GraphPath path(graph_ptr, 3, {0, 1, 2}, 10);
  ASSERT_FALSE(path.isValid());
}

TEST_F(DeletionGraph, PathWithUnorderedNodesIsInvalid) {
  GraphPath path(graph_ptr, 3, {2, 1}, 1);
  ASSERT_FALSE(path.isValid());
}

TEST_F(DeletionGraph, SingleNodePathWithEndBeforeStartIsInvalid) {
  GraphPath path(graph_ptr, 3, {0}, 1);
  ASSERT_FALSE(path.isValid());
}

TEST_F(SwapGraph, DisconnectedPathIsInvalid) {
  GraphPath path(graph_ptr, 0, {0, 3}, 0);
  ASSERT_FALSE(path.isValid());
}

TEST_F(DeletionGraph, PathIsEncodedAsString) {
  {
    GraphPath path(graph_ptr, 0, {0}, 1);
    ASSERT_EQ("(0@0)-(0@1)", path.encode());
  }
  {
    GraphPath path(graph_ptr, 3, {0, 1, 2}, 1);
    ASSERT_EQ("(0@3)-(1)-(2@1)", path.encode());
  }
}

TEST_F(DeletionGraph, PathStartPositionIsExtended) {
  {
    GraphPath path(graph_ptr, 3, {0, 1}, 1);
    GraphPath extended_path = path.extendStartPosition(3);
    GraphPath expected_path(graph_ptr, 0, {0, 1}, 1);
    ASSERT_EQ(expected_path, extended_path);
  }
  {
    GraphPath path(graph_ptr, 3, {0, 1}, 1);
    GraphPath extended_path = path.extendEndPosition(2);
    GraphPath expected_path(graph_ptr, 3, {0, 1}, 3);
    ASSERT_EQ(expected_path, extended_path);
  }
}

TEST_F(DeletionGraph, InvalidPathExtensionCausesError) {
  GraphPath path(graph_ptr, 3, {0, 1}, 1);
  ASSERT_ANY_THROW(path.extendStartPosition(4));
  ASSERT_ANY_THROW(path.extendEndPosition(4));
}

TEST_F(SwapGraph, PathIsExtendedToAdjacentNode) {
  {
    GraphPath path(graph_ptr, 1, {1, 3}, 2);
    GraphPath extended_path = path.extendStartNodeTo(0);
    const int32_t new_left_pos = left_flank.length() - 1;
    GraphPath expected_path(graph_ptr, new_left_pos, {0, 1, 3}, 2);
    ASSERT_EQ(expected_path, extended_path);
  }
  {
    GraphPath path(graph_ptr, 1, {0}, 2);
    GraphPath extended_path = path.extendEndNodeTo(1);
    GraphPath expected_path(graph_ptr, 1, {0, 1}, 0);
    ASSERT_EQ(expected_path, extended_path);
  }
}

TEST_F(SwapGraph, ExtendingPathToNonadjacentNodeCausesError) {
  {
    GraphPath path(graph_ptr, 1, {2, 3}, 3);
    ASSERT_ANY_THROW(path.extendStartNodeTo(1));
  }
  {
    GraphPath path(graph_ptr, 1, {0}, 3);
    ASSERT_ANY_THROW(path.extendEndNodeTo(3));
  }
}

TEST_F(DeletionGraph, EndOfPathIsExtendedBySpecifiedLength) {
  GraphPath path(graph_ptr, 2, {0}, 4);
  {
    const int32_t start_extension = 0;
    const int32_t end_extension = 1;
    const list<GraphPath> path_extensions =
        path.extendBy(start_extension, end_extension);
    const list<GraphPath> expected_path_extensions = {
        GraphPath(graph_ptr, 2, {0}, 5)};
    ASSERT_EQ(expected_path_extensions, path_extensions);
  }

  {
    const int32_t start_extension = 0;
    const int32_t end_extension = 2;
    const list<GraphPath> path_extensions =
        path.extendBy(start_extension, end_extension);
    const list<GraphPath> expected_path_extensions = {
        GraphPath(graph_ptr, 2, {0, 1}, 0), GraphPath(graph_ptr, 2, {0, 2}, 0)};
    ASSERT_EQ(expected_path_extensions, path_extensions);
  }
}

TEST_F(SwapGraph, EndOfPathIsExtendedBySpecifiedLength) {
  const GraphPath path(graph_ptr, 0, {0}, 4);
  const int32_t start_extension = 0;
  const int32_t end_extension = 5;
  const list<GraphPath> path_extensions =
      path.extendBy(start_extension, end_extension);
  const list<GraphPath> expected_path_extensions = {
      GraphPath(graph_ptr, 0, {0, 1, 3}, 0),
      GraphPath(graph_ptr, 0, {0, 2, 3}, 1)};
  ASSERT_EQ(expected_path_extensions, path_extensions);
}

TEST_F(DeletionGraph, PathIsExtendedBySpecifiedLength) {
  GraphPath path(graph_ptr, 2, {0}, 4);
  {
    const int32_t start_extension = 2;
    const int32_t end_extension = 1;
    const list<GraphPath> path_extensions =
        path.extendBy(start_extension, end_extension);
    const list<GraphPath> expected_path_extensions = {
        GraphPath(graph_ptr, 0, {0}, 5)};
    ASSERT_EQ(expected_path_extensions, path_extensions);
  }
}

TEST_F(DoubleSwapGraph, PathIsExtendedBySpecifiedLength) {
  GraphPath path(graph_ptr, 1, {3}, 3);
  {
    const int32_t start_extension = 5;
    const int32_t end_extension = 1;
    const list<GraphPath> path_extensions =
        path.extendBy(start_extension, end_extension);
    const list<GraphPath> expected_path_extensions = {
        GraphPath(graph_ptr, 0, {1, 3, 4}, 0),
        GraphPath(graph_ptr, 0, {1, 3, 5}, 0),
        GraphPath(graph_ptr, 4, {0, 2, 3, 4}, 0),
        GraphPath(graph_ptr, 4, {0, 2, 3, 5}, 0)};
    ASSERT_EQ(expected_path_extensions, path_extensions);
  }
}

TEST_F(DoubleSwapGraph, PathExtensionsThatAreTooLongAreNotOutput) {
  GraphPath path(graph_ptr, 1, {3}, 3);
  {
    const int32_t start_extension = 50;
    const int32_t end_extension = 1;
    const list<GraphPath> path_extensions =
        path.extendBy(start_extension, end_extension);
    ASSERT_TRUE(path_extensions.empty());
  }
  {
    const int32_t start_extension = 10;
    const int32_t end_extension = 9;
    const list<GraphPath> path_extensions =
        path.extendBy(start_extension, end_extension);
    const list<GraphPath> expected_path_extensions = {
        GraphPath(graph_ptr, 0, {0, 1, 3, 4, 6}, 3)};
    ASSERT_EQ(expected_path_extensions, path_extensions);
  }
}*/
