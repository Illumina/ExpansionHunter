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

TEST(TraversingPath, TypicalPath_NodeIdsTraversed) {
  GraphSharedPtr graph_ptr = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");
  GraphPath path(graph_ptr, 3, {1, 2}, 1);

  vector<int32_t> node_ids;
  for (int32_t node_id : path) {
    node_ids.push_back(node_id);
  }

  ASSERT_EQ(path.NodeIds(), node_ids);
}

TEST(GettingPathSequence, TypicalPathOnDeletionGraph_SequenceReturned) {
  GraphSharedPtr graph_ptr = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");

  {
    GraphPath path(graph_ptr, 3, {0}, 3);
    EXPECT_EQ("A", path.Seq());
  }

  {
    GraphPath path(graph_ptr, 3, {1}, 4);
    EXPECT_EQ("GG", path.Seq());
  }

  {
    GraphPath path(graph_ptr, 3, {0, 1, 2}, 1);
    EXPECT_EQ("ACCTTTGGAT", path.Seq());
  }
}

TEST(GettingPathSequence, TypicalPathOnStrGraph_SequenceReturned) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");
  GraphPath path(graph_ptr, 1, {0, 1, 1, 2}, 0);
  EXPECT_EQ("TTATATC", path.Seq());
}

TEST(CheckingIfPathOverlapsNode, TypicalPath_OverlapChecked) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");
  GraphPath path(graph_ptr, 1, {1, 1, 2}, 0);
  EXPECT_TRUE(path.OverlapsNode(1));
  EXPECT_FALSE(path.OverlapsNode(0));
}

TEST(GettingLengthOfPathOverEachNode, TypicalPathOnStrGraph_LengthReturned) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  GraphPath path(graph_ptr, 2, {0, 1, 1}, 0);

  EXPECT_EQ((size_t)1, path.GetOverlapWithNodeByIndex(0));
  EXPECT_EQ((size_t)2, path.GetOverlapWithNodeByIndex(1));
  EXPECT_EQ((size_t)1, path.GetOverlapWithNodeByIndex(2));
}

TEST(GettingLengthOfPathOverEachNode, IndexOutOfBounds_ExceptionRaised) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  GraphPath path(graph_ptr, 2, {0, 1, 1}, 0);

  EXPECT_ANY_THROW(path.GetOverlapWithNodeByIndex(-1));
  EXPECT_ANY_THROW(path.GetOverlapWithNodeByIndex(3));
}

TEST(GettingPathLength, TypicalPathOnStrGraph_LengthReturned) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  {
    GraphPath path(graph_ptr, 2, {0}, 2);
    EXPECT_EQ((size_t)1, path.Length());
  }

  {
    GraphPath path(graph_ptr, 0, {1}, 1);
    EXPECT_EQ((size_t)2, path.Length());
  }

  {
    GraphPath path(graph_ptr, 2, {0, 1, 1}, 0);
    EXPECT_EQ((size_t)4, path.Length());
  }
}

TEST(GettingPathSequenceOnNode, TypicalPathOnStrGraph_SequenceReturned) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  {
    GraphPath path(graph_ptr, 1, {0, 1, 1, 2}, 0);
    EXPECT_EQ("TT", path.SeqOnNodeByIndex(0));
    EXPECT_EQ("AT", path.SeqOnNodeByIndex(1));
    EXPECT_EQ("AT", path.SeqOnNodeByIndex(2));
    EXPECT_EQ("C", path.SeqOnNodeByIndex(3));
  }

  {
    GraphPath path(graph_ptr, 1, {1, 1}, 1);
    EXPECT_EQ("T", path.SeqOnNodeByIndex(0));
    EXPECT_EQ("AT", path.SeqOnNodeByIndex(1));
  }
}

TEST(ValidatingPath, WellFormedPath_IsValid) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  GraphPath path(graph_ptr, 1, {0, 1, 1, 2}, 0);
  ASSERT_TRUE(path.IsValid());
}

TEST(ValidatingPath, PathStartingOutsideOfNodeSequence_IsInvalid) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  GraphPath path(graph_ptr, 3, {0, 1, 2}, 1);
  ASSERT_FALSE(path.IsValid());
}

TEST(ValidatingPath, PathEndingOutsideOfNodeSequence_IsInvalid) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  GraphPath path(graph_ptr, 3, {0, 1, 2}, 10);
  ASSERT_FALSE(path.IsValid());
}

TEST(ValidatingPath, PathWithUnorderedNodes_IsInvalid) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  GraphPath path(graph_ptr, 1, {2, 1}, 1);
  ASSERT_FALSE(path.IsValid());
}

TEST(ValidatingPath, SingleNodePathWithEndBeforeStart_IsInvalid) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  GraphPath path(graph_ptr, 3, {0}, 1);
  ASSERT_FALSE(path.IsValid());
}

TEST(ValidatingPath, DisconnectedPath_IsInvalid) {
  GraphSharedPtr graph_ptr = MakeSwapGraph("TTT", "AT", "GG", "CCCCC");

  GraphPath path(graph_ptr, 0, {0, 3}, 0);
  ASSERT_FALSE(path.IsValid());
}

TEST(EncodingPaths, TypicalPath_EncodedAsString) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  {
    GraphPath path(graph_ptr, 0, {0}, 1);
    ASSERT_EQ("(0@0)-(0@1)", path.Encode());
  }

  {
    GraphPath path(graph_ptr, 1, {0, 1, 1, 2}, 0);
    ASSERT_EQ("(0@1)-(1)-(1)-(2@0)", path.Encode());
  }
}

TEST(MovePathAlongNode, TypicalPath_StartPositionMoved) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");
  GraphPath shorter_path(graph_ptr, 2, {0, 1}, 1);
  GraphPath longer_path(graph_ptr, 0, {0, 1}, 1);

  EXPECT_EQ(longer_path, shorter_path.MoveStartBy(2));
  EXPECT_EQ(shorter_path, longer_path.MoveStartBy(-2));
}

TEST(MovePathAlongNode, TypicalPath_EndPositionMoved) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");
  GraphPath shorter_path(graph_ptr, 1, {0, 1, 1}, 0);
  GraphPath longer_path(graph_ptr, 1, {0, 1, 1}, 1);

  EXPECT_EQ(longer_path, shorter_path.MoveEndBy(1));
  EXPECT_EQ(shorter_path, longer_path.MoveEndBy(-1));
}

TEST(MovePathAlongNode, ExtensionPastNodeBoundaries_ExceptionRaised) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  GraphPath path(graph_ptr, 2, {0, 1}, 1);
  EXPECT_ANY_THROW(path.MoveStartBy(3));
  EXPECT_ANY_THROW(path.MoveStartBy(-1));
  EXPECT_ANY_THROW(path.MoveEndBy(1));
  EXPECT_ANY_THROW(path.MoveEndBy(-2));
}

TEST(ExtendingPathToNode, TypicalPathInSwapGraph_StartPositionMoved) {
  GraphSharedPtr graph_ptr = MakeSwapGraph("TTT", "AT", "GG", "CCCCC");

  GraphPath path(graph_ptr, 1, {1, 3}, 2);

  {
    GraphPath expected_path(graph_ptr, 2, {0, 1, 3}, 2);
    EXPECT_EQ(expected_path, path.ExtendStartToNode(0));
  }

  {
    GraphPath expected_path(graph_ptr, 0, {3}, 2);
    EXPECT_EQ(expected_path, path.RemoveStartNode());
  }
}

TEST(ExtendingPathToNode, TypicalPathInSwapGraph_EndPositionMoved) {
  GraphSharedPtr graph_ptr = MakeSwapGraph("TTT", "AT", "GG", "CCCCC");

  GraphPath path(graph_ptr, 1, {0, 2}, 1);
  {
    GraphPath expected_path(graph_ptr, 1, {0, 2, 3}, 0);
    ASSERT_EQ(expected_path, path.ExtendEndToNode(3));
  }

  {
    GraphPath expected_path(graph_ptr, 1, {0}, 2);
    ASSERT_EQ(expected_path, path.RemoveEndNode());
  }
}

TEST(ExtendingPathToNode, ExtendingPathToNonadjacentNode_ExceptionThrown) {
  GraphSharedPtr graph_ptr = MakeSwapGraph("TTT", "AT", "GG", "CCCCC");

  {
    GraphPath path(graph_ptr, 1, {2, 3}, 1);
    EXPECT_ANY_THROW(path.ExtendStartToNode(1));
  }

  {
    GraphPath path(graph_ptr, 1, {0}, 2);
    EXPECT_ANY_THROW(path.ExtendEndToNode(3));
  }
}

TEST(ExtendingPathsByGivenLength, TypicalPathInStrGraph_PathExtended) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AT", "CCCCC");

  {
    GraphPath path(graph_ptr, 1, {0}, 1);
    const int32_t start_extension = 0;
    const int32_t end_extension = 6;
    const list<GraphPath> path_extensions =
        path.ExtendBy(start_extension, end_extension);

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
        path.ExtendBy(start_extension, end_extension);

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
  GraphSharedPtr graph_ptr = MakeStrGraph("T", "A", "C");

  GraphPath path(graph_ptr, 0, {1}, 0);
  const int32_t start_extension = 3;
  const int32_t end_extension = 3;
  const list<GraphPath> path_extensions =
      path.ExtendBy(start_extension, end_extension);

  const list<GraphPath> expected_path_extensions = {
      GraphPath(graph_ptr, 0, {0, 1, 1, 1, 1, 1, 1}, 0),
      GraphPath(graph_ptr, 0, {0, 1, 1, 1, 1, 1, 2}, 0),
      GraphPath(graph_ptr, 0, {1, 1, 1, 1, 1, 1, 1}, 0),
      GraphPath(graph_ptr, 0, {1, 1, 1, 1, 1, 1, 2}, 0)};
  ASSERT_EQ(expected_path_extensions, path_extensions);
}

TEST(ShrinkingPathEnds, TypicalPathInStrGraph_EndShrank) {
  GraphSharedPtr graph_ptr = MakeStrGraph("ATAT", "C", "CCTT");
  GraphPath path(graph_ptr, 2, {0, 1, 2}, 2);

  GraphPath shrank_path = path.ShrinkEndBy(4);

  GraphPath expected_path(graph_ptr, 2, {0}, 3);
  ASSERT_EQ(expected_path, shrank_path);
}

TEST(ShrinkingPathsByGivenLength, TypicalPathInStrGraph_PathShrank) {
  GraphSharedPtr graph_ptr = MakeStrGraph("TTT", "AC", "CCC");

  GraphPath path(graph_ptr, 1, {0, 1, 1, 2}, 2);
  const int32_t start_shrink_len = 5;
  const int32_t end_shrink_len = 3;
  const GraphPath shrank_path = path.ShrinkBy(start_shrink_len, end_shrink_len);

  const GraphPath expected_path(graph_ptr, 1, {1}, 1);
  ASSERT_EQ(expected_path, shrank_path);
}