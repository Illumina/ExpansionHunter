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

#include "graphcore/Path.hh"

#include "gtest/gtest.h"

#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"

using std::string;
using std::vector;

using namespace graphtools;

TEST(CreatingPath, WellFormedPath_NoExceptionThrown)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    ASSERT_NO_THROW(Path(&graph, 1, { 0, 1, 1, 2 }, 0));
}

TEST(CreatingPath, ZeroLengthPathSpanningAnEdge_NoExceptionThrown)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    ASSERT_NO_THROW(Path(&graph, 3, { 0, 1, 1, 2 }, 0));
}

#ifdef _DEBUG

TEST(CreatingPath, PathWithUnorderedNodes_ExceptionThrown)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    ASSERT_ANY_THROW(Path(&graph, 1, { 2, 1 }, 1));
}

TEST(CreatingPath, PathStartingOutsideOfNodesequence_ExceptionThrown)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    ASSERT_ANY_THROW(Path(&graph, 4, { 0, 1, 2 }, 1));
}

TEST(CreatingPath, PathEndingOutsideOfNodesequence_ExceptionThrown)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    ASSERT_ANY_THROW(Path(&graph, 3, { 0, 1, 2 }, 10));
}

TEST(CreatingPath, PathWithEndBeforeStart_ExceptionThrown)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    ASSERT_ANY_THROW(Path(&graph, 3, { 0 }, 1));
}

TEST(CreatingPath, DisconnectedPath_ExceptionThrown)
{
    Graph graph = makeSwapGraph("TTT", "AT", "GG", "CCCCC");
    ASSERT_ANY_THROW(Path(&graph, 0, { 0, 3 }, 0));
}

#endif

TEST(TraversingPath, TypicalPath_NodeIdsTraversed)
{
    Graph graph = makeDeletionGraph("AAAACC", "TTTGG", "ATTT");
    Path path(&graph, 3, { 1, 2 }, 1);

    vector<NodeId> node_ids;
    for (NodeId node_id : path)
    {
        node_ids.push_back(node_id);
    }

    ASSERT_EQ(path.nodeIds(), node_ids);
}

TEST(GettingPathsequence, TypicalPathOnDeletionGraph_SequenceReturned)
{
    Graph graph = makeDeletionGraph("AAAACC", "TTTGG", "ATTT");

    {
        Path path(&graph, 3, { 0 }, 3);
        EXPECT_EQ("", path.seq());
    }

    {
        Path path(&graph, 3, { 1 }, 4);
        EXPECT_EQ("G", path.seq());
    }

    {
        Path path(&graph, 3, { 0, 1, 2 }, 1);
        EXPECT_EQ("ACCTTTGGA", path.seq());
    }
}

TEST(GettingPathsequence, TypicalPathOnStrGraph_sequenceReturned)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    Path path(&graph, 1, { 0, 1, 1, 2 }, 0);
    EXPECT_EQ("TTATAT", path.seq());
}

TEST(CheckingIfPathOverlapsNode, TypicalPath_OverlapChecked)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    Path path(&graph, 1, { 1, 1, 2 }, 0);
    EXPECT_TRUE(path.checkOverlapWithNode(1));
    EXPECT_TRUE(path.checkOverlapWithNode(2));
    EXPECT_FALSE(path.checkOverlapWithNode(0));
}

TEST(GettingPathBoundsOnNodeByIndex, TypicalPath_BoundsComputed)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

    {
        Path path(&graph, 1, { 0 }, 2);
        EXPECT_EQ(1, path.getStartPositionOnNodeByIndex(0));
        EXPECT_EQ(2, path.getEndPositionOnNodeByIndex(0));
    }

    {
        Path path(&graph, 1, { 1, 1, 2 }, 3);
        EXPECT_EQ(1, path.getStartPositionOnNodeByIndex(0));
        EXPECT_EQ(2, path.getEndPositionOnNodeByIndex(0));

        EXPECT_EQ(0, path.getStartPositionOnNodeByIndex(1));
        EXPECT_EQ(2, path.getEndPositionOnNodeByIndex(1));

        EXPECT_EQ(0, path.getStartPositionOnNodeByIndex(2));
        EXPECT_EQ(3, path.getEndPositionOnNodeByIndex(2));

        EXPECT_ANY_THROW(path.getStartPositionOnNodeByIndex(-1));
        EXPECT_ANY_THROW(path.getEndPositionOnNodeByIndex(-1));

        EXPECT_ANY_THROW(path.getStartPositionOnNodeByIndex(3));
        EXPECT_ANY_THROW(path.getEndPositionOnNodeByIndex(3));
    }
}

TEST(GettingLengthOfPathOverEachNode, TypicalPathOnStrGraph_LengthReturned)
{
    {
        Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

        Path path(&graph, 2, { 0, 1, 1 }, 0);

        EXPECT_EQ(1ul, path.getNodeOverlapLengthByIndex(0));
        EXPECT_EQ(2ul, path.getNodeOverlapLengthByIndex(1));
        EXPECT_EQ(0ul, path.getNodeOverlapLengthByIndex(2));
    }

    {
        Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

        Path path(&graph, 3, { 0, 1, 1, 2 }, 5);

        EXPECT_EQ(0ul, path.getNodeOverlapLengthByIndex(0));
        EXPECT_EQ(2ul, path.getNodeOverlapLengthByIndex(1));
        EXPECT_EQ(2ul, path.getNodeOverlapLengthByIndex(2));
        EXPECT_EQ(5ul, path.getNodeOverlapLengthByIndex(3));
    }
}

TEST(GettingLengthOfPathOverEachNode, IndexOutOfBounds_ExceptionRaised)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

    Path path(&graph, 2, { 0, 1, 1 }, 0);

    EXPECT_ANY_THROW(path.getNodeOverlapLengthByIndex(-1));
    EXPECT_ANY_THROW(path.getNodeOverlapLengthByIndex(3));
}

TEST(GettingPathLength, TypicalPathOnStrGraph_LengthReturned)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

    {
        Path path(&graph, 2, { 0 }, 2);
        EXPECT_EQ(0ul, path.length());
    }

    {
        Path path(&graph, 0, { 1 }, 1);
        EXPECT_EQ(1ul, path.length());
    }

    {
        Path path(&graph, 2, { 0, 1, 1 }, 0);
        EXPECT_EQ(3ul, path.length());
    }

    {
        Path path(&graph, 3, { 0, 1, 1 }, 0);
        EXPECT_EQ(2ul, path.length());
    }
}

TEST(GettingPathsequenceOnNode, TypicalPathOnStrGraph_sequenceReturned)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

    {
        Path path(&graph, 1, { 0, 1, 1, 2 }, 0);
        EXPECT_EQ("TT", path.getNodeSeq(0));
        EXPECT_EQ("AT", path.getNodeSeq(1));
        EXPECT_EQ("AT", path.getNodeSeq(2));
        EXPECT_EQ("", path.getNodeSeq(3));
    }

    {
        Path path(&graph, 1, { 1, 1 }, 1);
        EXPECT_EQ("T", path.getNodeSeq(0));
        EXPECT_EQ("A", path.getNodeSeq(1));
    }
}

TEST(EncodingPaths, TypicalPath_EncodedAsString)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

    {
        Path path(&graph, 0, { 0 }, 1);
        ASSERT_EQ("(0@0)-(0@1)", path.encode());
    }

    {
        Path path(&graph, 1, { 0, 1, 1, 2 }, 0);
        ASSERT_EQ("(0@1)-(1)-(1)-(2@0)", path.encode());
    }
}

TEST(MovePathAlongNode, TypicalPath_StartPositionMoved)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    Path shorter_path(&graph, 3, { 0, 1 }, 1);
    Path longer_path(&graph, 0, { 0, 1 }, 1);

    {
        Path path(shorter_path);
        path.shiftStartAlongNode(3);
        EXPECT_EQ(longer_path, path);
    }

    {
        Path path(longer_path);
        path.shiftStartAlongNode(-3);
        EXPECT_EQ(shorter_path, path);
    }
}

TEST(MovePathAlongNode, TypicalPath_EndPositionMoved)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    Path shorter_path(&graph, 1, { 0, 1, 1 }, 0);
    Path longer_path(&graph, 1, { 0, 1, 1 }, 1);

    {
        Path path(shorter_path);
        path.shiftEndAlongNode(1);
        EXPECT_EQ(longer_path, path);
    }

    {
        Path path(longer_path);
        path.shiftEndAlongNode(-1);
        EXPECT_EQ(shorter_path, path);
    }
}

TEST(MovePathAlongNode, ExtensionPastNodeBoundaries_ExceptionRaised)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

    {
        Path path(&graph, 2, { 0, 1 }, 1);
        EXPECT_ANY_THROW(path.shiftStartAlongNode(3));
    }

    {
        Path path(&graph, 2, { 0, 1 }, 1);
        EXPECT_ANY_THROW(path.shiftStartAlongNode(-2));
    }

    {
        Path path(&graph, 2, { 0, 1 }, 1);
        EXPECT_ANY_THROW(path.shiftEndAlongNode(2));
    }

    {
        Path path(&graph, 2, { 0, 1 }, 1);
        EXPECT_ANY_THROW(path.shiftEndAlongNode(-2));
    }
}

TEST(ExtendingPathToNode, TypicalPathInSwapGraph_StartPositionMoved)
{
    Graph graph = makeSwapGraph("TTT", "AT", "GG", "CCCCC");

    {
        Path path(&graph, 1, { 1, 3 }, 2);
        path.extendStartToNode(0);
        Path expected_path(&graph, 3, { 0, 1, 3 }, 2);
        EXPECT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 1, { 1, 3 }, 2);
        path.removeStartNode();
        Path expected_path(&graph, 0, { 3 }, 2);
        EXPECT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 1, { 1, 3 }, 2);
        path.extendStartToIncludeNode(0);
        Path expected_path(&graph, 0, { 0, 1, 3 }, 2);
        EXPECT_EQ(expected_path, path);
    }
}

TEST(ExtendingPathToNode, TypicalPathInSwapGraph_EndPositionMoved)
{
    Graph graph = makeSwapGraph("TTT", "AT", "GG", "CCCCC");

    {
        Path path(&graph, 1, { 0, 2 }, 1);
        path.extendEndToNode(3);
        Path expected_path(&graph, 1, { 0, 2, 3 }, 0);
        EXPECT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 1, { 0, 2 }, 1);
        path.removeEndNode();
        Path expected_path(&graph, 1, { 0 }, 3);
        EXPECT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 1, { 0, 2 }, 1);
        path.extendEndToIncludeNode(3);
        Path expected_path(&graph, 1, { 0, 2, 3 }, 5);
        EXPECT_EQ(expected_path, path);
    }
}

TEST(ExtendingPathToNode, ExtendingPathToNonadjacentNode_ExceptionThrown)
{
    Graph graph = makeSwapGraph("TTT", "AT", "GG", "CCCCC");

    {
        Path path(&graph, 1, { 2, 3 }, 1);
        EXPECT_ANY_THROW(path.extendStartToNode(1));
    }

    {
        Path path(&graph, 1, { 0 }, 2);
        EXPECT_ANY_THROW(path.extendEndToNode(3));
    }
}

TEST(RemovingZeroLengthStarts, TypicalPaths_StartRemovedIfAppropriate)
{
    Graph graph = makeStrGraph("ATAT", "C", "CCTT");

    {
        Path path(&graph, 4, { 0, 1, 2 }, 2);
        path.removeZeroLengthStart();

        Path expected_path(&graph, 0, { 1, 2 }, 2);

        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 3, { 0, 1, 2 }, 2);
        path.removeZeroLengthStart();

        Path expected_path(&graph, 3, { 0, 1, 2 }, 2);

        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 4, { 0 }, 4);
        path.removeZeroLengthStart();

        Path expected_path(&graph, 4, { 0 }, 4);

        ASSERT_EQ(expected_path, path);
    }
}

TEST(RemovingZeroLengthEnds, TypicalPaths_EndRemovedIfAppropriate)
{
    Graph graph = makeStrGraph("ATAT", "C", "CCTT");

    {
        Path path(&graph, 0, { 0, 1, 2 }, 0);
        path.removeZeroLengthEnd();

        Path expected_path(&graph, 0, { 0, 1 }, 1);

        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 0, { 0, 1, 2 }, 1);
        path.removeZeroLengthEnd();

        Path expected_path(&graph, 0, { 0, 1, 2 }, 1);

        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 4, { 0 }, 4);
        path.removeZeroLengthEnd();

        Path expected_path(&graph, 4, { 0 }, 4);

        ASSERT_EQ(expected_path, path);
    }
}

TEST(ShrinkingStartOfPath, TypicalPathInStrGraph_StartShrank)
{
    Graph graph = makeStrGraph("ATAT", "C", "CCTT");

    {
        Path path(&graph, 2, { 0, 1, 2 }, 2);
        path.shrinkStartBy(2);

        Path expected_path(&graph, 0, { 1, 2 }, 2);
        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 2, { 0, 1 }, 1);
        path.shrinkStartBy(3);

        Path expected_path(&graph, 1, { 1 }, 1);
        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 4, { 0, 1 }, 1);
        path.shrinkStartBy(1);

        Path expected_path(&graph, 1, { 1 }, 1);
        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 4, { 0, 1 }, 1);
        path.shrinkStartBy(0);

        Path expected_path(&graph, 0, { 1 }, 1);
        ASSERT_EQ(expected_path, path);
    }
}

TEST(ShrinkingEndOfPath, TypicalPathInStrGraph_EndShrank)
{
    Graph graph = makeStrGraph("ATAT", "C", "CCTT");

    {
        Path path(&graph, 2, { 0, 1, 2 }, 2);
        path.shrinkEndBy(3);

        Path expected_path(&graph, 2, { 0 }, 4);
        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 0, { 1, 2 }, 2);
        path.shrinkEndBy(3);

        Path expected_path(&graph, 0, { 1 }, 0);
        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 0, { 1, 2 }, 0);
        path.shrinkEndBy(1);

        Path expected_path(&graph, 0, { 1 }, 0);
        ASSERT_EQ(expected_path, path);
    }

    {
        Path path(&graph, 0, { 1, 2 }, 0);
        path.shrinkEndBy(0);

        Path expected_path(&graph, 0, { 1 }, 1);
        ASSERT_EQ(expected_path, path);
    }
}

TEST(ShrinkingPathEnds, PathWithLoop_PathShrank)
{
    Graph graph = makeStrGraph("ATA", "CG", "TATTTTTTTTT");

    Path path(&graph, 1, { 0, 1, 1, 1, 2 }, 3);
    path.shrinkEndBy(5);

    Path expected_path(&graph, 1, { 0, 1, 1 }, 2);
    ASSERT_EQ(expected_path, path);
}

TEST(ShrinkingPathsByGivenLength, TypicalPathInStrGraph_PathShrank)
{
    Graph graph = makeStrGraph("TTT", "AC", "CCC");

    Path path(&graph, 1, { 0, 1, 1, 2 }, 2);
    const int32_t start_shrink_len = 5;
    const int32_t end_shrink_len = 3;
    path.shrinkBy(start_shrink_len, end_shrink_len);

    const Path expected_path(&graph, 1, { 1 }, 1);
    ASSERT_EQ(expected_path, path);
}

TEST(ComputingPathDistance, DistanceFromStart_DistanceReturned)
{
    Graph graph = makeDeletionGraph("TTT", "AC", "CCC");

    Path path(&graph, 1, { 0, 1, 2 }, 2);
    ASSERT_EQ(0, path.getDistanceFromPathStart(0, 1));
    ASSERT_EQ(4, path.getDistanceFromPathStart(2, 0));
    ASSERT_EQ(3, path.getDistanceFromPathStart(1, 1));
}

TEST(ComputingPathDistance, DistanceFromStart_ExceptionWhenNodeNotInPath)
{
    Graph graph = makeDeletionGraph("TTT", "AC", "CCCC");

    Path path(&graph, 1, { 0, 2 }, 2);

    ASSERT_THROW(path.getDistanceFromPathStart(0, 0), std::logic_error);
    ASSERT_THROW(path.getDistanceFromPathStart(1, 0), std::logic_error);
    ASSERT_THROW(path.getDistanceFromPathStart(2, 3), std::logic_error);
}

TEST(ComparingPaths, TypicalPaths_Compared)
{
    Graph graph = makeStrGraph("TTT", "AC", "CCC");

    {
        Path path_a(&graph, 1, { 0, 1, 2 }, 1);
        Path path_b(&graph, 1, { 0, 1, 2 }, 2);
        EXPECT_TRUE(path_a < path_b);
        EXPECT_FALSE(path_b < path_a);
        EXPECT_FALSE(path_a == path_b);
    }

    {
        Path path_a(&graph, 0, { 0, 1, 1 }, 1);
        Path path_b(&graph, 0, { 0, 1, 2 }, 1);
        EXPECT_TRUE(path_a < path_b);
        EXPECT_FALSE(path_b < path_a);
        EXPECT_FALSE(path_a == path_b);
    }
}
