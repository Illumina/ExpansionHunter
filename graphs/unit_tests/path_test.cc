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

TEST(GettingLengthOfPathOverEachNode, IndexOutOfBounds_ExceptionRaised) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 2, {0, 1, 1}, 0);

  EXPECT_ANY_THROW(path.GetOverlapWithNodeByIndex(-1));
  EXPECT_ANY_THROW(path.GetOverlapWithNodeByIndex(3));
}

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
}

TEST(GettingPathSequenceOnNode, TypicalPathOnStrGraph_SequenceReturned) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 1, {0, 1, 1, 2}, 0);
  EXPECT_EQ("TT", path.SeqOnNodeByIndex(0));
  EXPECT_EQ("AT", path.SeqOnNodeByIndex(1));
  EXPECT_EQ("AT", path.SeqOnNodeByIndex(2));
  EXPECT_EQ("C", path.SeqOnNodeByIndex(3));
}

TEST(ValidatingPath, WellFormedPath_IsValid) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 1, {0, 1, 1, 2}, 0);
  ASSERT_TRUE(path.isValid());
}

TEST(ValidatingPath, PathStartingOutsideOfNodeSequence_IsInvalid) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 3, {0, 1, 2}, 1);
  ASSERT_FALSE(path.isValid());
}

TEST(ValidatingPath, PathEndingOutsideOfNodeSequence_IsInvalid) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 3, {0, 1, 2}, 10);
  ASSERT_FALSE(path.isValid());
}

TEST(ValidatingPath, PathWithUnorderedNodes_IsInvalid) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 1, {2, 1}, 1);
  ASSERT_FALSE(path.isValid());
}

TEST(ValidatingPath, SingleNodePathWithEndBeforeStart_IsInvalid) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 3, {0}, 1);
  ASSERT_FALSE(path.isValid());
}

TEST(ValidatingPath, DisconnectedPath_IsInvalid) {
  Graph graph = makeSwapGraph("TTT", "AT", "GG", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 0, {0, 3}, 0);
  ASSERT_FALSE(path.isValid());
}

TEST(EncodingPaths, TypicalPath_EncodedAsString) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  {
    GraphPath path(graph_ptr, 0, {0}, 1);
    ASSERT_EQ("(0@0)-(0@1)", path.encode());
  }

  {
    GraphPath path(graph_ptr, 1, {0, 1, 1, 2}, 0);
    ASSERT_EQ("(0@1)-(1)-(1)-(2@0)", path.encode());
  }
}

TEST(ExtendingPathAlongNode, TypicalPath_StartPositionExtended) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  {
    GraphPath path(graph_ptr, 2, {0, 1}, 1);
    GraphPath extended_path = path.extendStartPosition(2);
    GraphPath expected_path(graph_ptr, 0, {0, 1}, 1);
    ASSERT_EQ(expected_path, extended_path);
  }

  {
    GraphPath path(graph_ptr, 1, {0, 1, 1}, 0);
    GraphPath extended_path = path.extendEndPosition(1);
    GraphPath expected_path(graph_ptr, 1, {0, 1, 1}, 1);
    ASSERT_EQ(expected_path, extended_path);
  }
}

TEST(ExtendingPathAlongNode, ExtensionPastNodeBoundaries_ExceptionRaised) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 2, {0, 1}, 1);
  EXPECT_ANY_THROW(path.extendStartPosition(3));
  EXPECT_ANY_THROW(path.extendEndPosition(1));
}

TEST(ExtendingPathToNode, TypicalPathInSwapGraph_PathExtended) {
  Graph graph = makeSwapGraph("TTT", "AT", "GG", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  {
    GraphPath path(graph_ptr, 1, {1, 3}, 2);
    GraphPath expected_path(graph_ptr, 2, {0, 1, 3}, 2);
    ASSERT_EQ(expected_path, path.extendStartNodeTo(0));
  }

  {
    GraphPath path(graph_ptr, 1, {0}, 2);
    GraphPath expected_path(graph_ptr, 1, {0, 1}, 0);
    ASSERT_EQ(expected_path, path.extendEndNodeTo(1));
  }
}

TEST(ExtendingPathToNode, ExtendingPathToNonadjacentNode_ExceptionThrown) {
  Graph graph = makeSwapGraph("TTT", "AT", "GG", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  {
    GraphPath path(graph_ptr, 1, {2, 3}, 1);
    EXPECT_ANY_THROW(path.extendStartNodeTo(1));
  }
  {
    GraphPath path(graph_ptr, 1, {0}, 2);
    EXPECT_ANY_THROW(path.extendEndNodeTo(3));
  }
}

TEST(ExtendingPathsByGivenLength, TypicalPathInStrGraph_PathExtended) {
  Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  {
    GraphPath path(graph_ptr, 1, {0}, 1);
    const int32_t start_extension = 0;
    const int32_t end_extension = 6;
    const list<GraphPath> path_extensions =
        path.extendBy(start_extension, end_extension);

    const list<GraphPath> expected_path_extensions = {
        GraphPath(graph_ptr, 1, {0, 1, 1, 1}, 0),
        GraphPath(graph_ptr, 1, {0, 1, 1, 2}, 0),
        GraphPath(graph_ptr, 1, {0, 1, 2}, 2),
        GraphPath(graph_ptr, 1, {0, 2}, 4)};
    ASSERT_EQ(expected_path_extensions, path_extensions);
  }

  {
    GraphPath path(graph_ptr, 0, {1}, 1);
    const int32_t start_extension = 1;
    const int32_t end_extension = 1;
    const list<GraphPath> path_extensions =
        path.extendBy(start_extension, end_extension);

    const list<GraphPath> expected_path_extensions = {
        GraphPath(graph_ptr, 2, {0, 1, 1}, 0),
        GraphPath(graph_ptr, 2, {0, 1, 2}, 0),
        GraphPath(graph_ptr, 1, {1, 1, 1}, 0),
        GraphPath(graph_ptr, 1, {1, 1, 2}, 0),
    };
    ASSERT_EQ(expected_path_extensions, path_extensions);
  }
}

TEST(ExtendingPathsByGivenLength, TypicalPathInHomopolymerGraph_PathExtended) {
  Graph graph = makeStrGraph("T", "A", "C");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);

  GraphPath path(graph_ptr, 0, {1}, 0);
  const int32_t start_extension = 3;
  const int32_t end_extension = 3;
  const list<GraphPath> path_extensions =
      path.extendBy(start_extension, end_extension);

  const list<GraphPath> expected_path_extensions = {
      GraphPath(graph_ptr, 0, {0, 1, 1, 1, 1, 1, 1}, 0),
      GraphPath(graph_ptr, 0, {0, 1, 1, 1, 1, 1, 2}, 0),
      GraphPath(graph_ptr, 0, {1, 1, 1, 1, 1, 1, 1}, 0),
      GraphPath(graph_ptr, 0, {1, 1, 1, 1, 1, 1, 2}, 0)};
  ASSERT_EQ(expected_path_extensions, path_extensions);
}
