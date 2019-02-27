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

#include "graphcore/PathOperations.hh"

#include "gtest/gtest.h"

#include <list>
#include <string>
#include <vector>

#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/Path.hh"

using std::list;
using std::string;
using std::vector;

using namespace graphtools;

TEST(ExtendingPathStarts, TypicalPath_StartExtended)
{
    Graph graph = makeDeletionGraph("AAACC", "TTGGG", "TTAAA");

    {
        Path path(&graph, 4, { 0 }, 4);
        list<Path> extended_paths = extendPathStart(path, 1);

        list<Path> expected_paths = { Path(&graph, 3, { 0 }, 4) };
        ASSERT_EQ(expected_paths, extended_paths);
    }

    {
        Path path(&graph, 5, { 0, 2 }, 0);
        list<Path> extended_paths = extendPathStart(path, 2);

        list<Path> expected_paths = { Path(&graph, 3, { 0, 2 }, 0) };
        ASSERT_EQ(expected_paths, extended_paths);
    }

    {
        Path path(&graph, 0, { 2 }, 0);
        list<Path> extended_paths = extendPathStart(path, 2);

        list<Path> expected_paths = { Path(&graph, 3, { 0, 2 }, 0), Path(&graph, 3, { 1, 2 }, 0) };
        ASSERT_EQ(expected_paths, extended_paths);
    }
}

TEST(ExtendingPathEnds, TypicalPath_EndExtended)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

    Path path(&graph, 0, { 0 }, 1);
    const list<Path> path_extensions = extendPathEnd(path, 6);

    const list<Path> expected_path_extensions
        = { Path(&graph, 0, { 0, 1, 1 }, 2), Path(&graph, 0, { 0, 1, 2 }, 2), Path(&graph, 0, { 0, 2 }, 4) };
    ASSERT_EQ(expected_path_extensions, path_extensions);
}

TEST(ExtendingPathsByGivenLength, TypicalPathInStrGraph_PathExtended)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

    Path path(&graph, 0, { 1 }, 2);
    const int32_t start_extension = 1;
    const int32_t end_extension = 1;
    const list<Path> path_extensions = extendPath(path, start_extension, end_extension);

    const list<Path> expected_path_extensions = { Path(&graph, 2, { 0, 1, 1 }, 1), Path(&graph, 2, { 0, 1, 2 }, 1),
                                                  Path(&graph, 1, { 1, 1, 1 }, 1), Path(&graph, 1, { 1, 1, 2 }, 1) };
    ASSERT_EQ(expected_path_extensions, path_extensions);
}

TEST(ExtendingPathsByGivenLength, TypicalPathInHomopolymerGraph_PathExtended)
{
    Graph graph = makeStrGraph("T", "A", "C");

    Path path(&graph, 0, { 1 }, 0);
    const int32_t start_extension = 3;
    const int32_t end_extension = 3;
    const list<Path> path_extensions = extendPath(path, start_extension, end_extension);

    const list<Path> expected_path_extensions
        = { Path(&graph, 0, { 0, 1, 1, 1, 1, 1 }, 1), Path(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 1),
            Path(&graph, 0, { 1, 1, 1, 1, 1, 1 }, 1), Path(&graph, 0, { 1, 1, 1, 1, 1, 2 }, 1) };
    ASSERT_EQ(expected_path_extensions, path_extensions);
}

TEST(ExtendingPathsMatching, TypicalPath_ExtendedWithinNode)
{
    Graph graph = makeDeletionGraph("AAACC", "TTGGG", "TTAAA");

    const Path path(&graph, 2, { 1 }, 2);
    const string query = "TTGGG";

    {
        size_t qpos = 2;
        const Path extended_path = extendPathStartMatching(path, query, qpos);
        const Path expected_path(&graph, 0, { 1 }, 2);
        ASSERT_EQ(expected_path, extended_path);
        ASSERT_EQ(0ull, qpos);
    }

    {
        const size_t qpos = 2;
        const Path extended_path = extendPathEndMatching(path, query, qpos);
        const Path expected_path(&graph, 2, { 1 }, 5);
        ASSERT_EQ(expected_path, extended_path);
        ASSERT_EQ(2ull, qpos);
    }

    {
        size_t qpos = 2;
        const Path extended_path = extendPathMatching(path, query, qpos);
        const Path expected_path(&graph, 0, { 1 }, 5);
        ASSERT_EQ(expected_path, extended_path);
        ASSERT_EQ(0ull, qpos);
    }
}

TEST(ExtendingPathsMatching, TypicalPath_ExtendedAcrossNodes)
{
    Graph graph = makeDeletionGraph("AAACC", "TTGGG", "TTAAA");

    const Path path(&graph, 2, { 1 }, 2);
    const string query = "CTTGGGT";

    {
        size_t qpos = 3;
        const Path extended_path = extendPathStartMatching(path, query, qpos);
        const Path expected_path(&graph, 4, { 0, 1 }, 2);
        ASSERT_EQ(expected_path, extended_path);
        ASSERT_EQ(0ull, qpos);
    }

    {
        const size_t qpos = 3;
        const Path extended_path = extendPathEndMatching(path, query, qpos);
        const Path expected_path(&graph, 2, { 1, 2 }, 1);
        ASSERT_EQ(expected_path, extended_path);
        ASSERT_EQ(3ull, qpos);
    }

    {
        size_t qpos = 3;
        const Path extended_path = extendPathMatching(path, query, qpos);
        const Path expected_path(&graph, 4, { 0, 1, 2 }, 1);
        ASSERT_EQ(expected_path, extended_path);
        ASSERT_EQ(0ull, qpos);
    }
}

TEST(ExtendingPathsMatching, TypicalPath_ExtendedWhenUniqMatch)
{
    {
        Graph graph = makeDeletionGraph("AAACC", "TTGGG", "TTAAA");
        const Path path(&graph, 4, { 0 }, 4);

        {
            const string query = "CTTGG";
            size_t qpos = 0;
            const Path extended_path = extendPathMatching(path, query, qpos);
            const Path expected_path(&graph, 4, { 0, 1 }, 4);
            ASSERT_EQ(expected_path, extended_path);
            ASSERT_EQ(0ull, qpos);
        }

        {
            const string query = "CTTAA";
            size_t qpos = 0;
            const Path extended_path = extendPathMatching(path, query, qpos);
            const Path expected_path(&graph, 4, { 0, 2 }, 4);
            ASSERT_EQ(expected_path, extended_path);
            ASSERT_EQ(qpos, 0ull);
        }
    }

    {
        Graph graph = makeDeletionGraph("AAACC", "ATGCC", "TTAAA");
        const Path path(&graph, 0, { 2 }, 0);

        {
            const string query = "TGCCT";
            size_t qpos = 4;
            const Path extended_path = extendPathMatching(path, query, qpos);
            const Path expected_path(&graph, 1, { 1, 2 }, 1);
            ASSERT_EQ(expected_path, extended_path);
            ASSERT_EQ(0ull, qpos);
        }

        {
            const string query = "AACCT";
            size_t qpos = 4;
            const Path extended_path = extendPathMatching(path, query, qpos);
            const Path expected_path(&graph, 1, { 0, 2 }, 1);
            ASSERT_EQ(expected_path, extended_path);
            ASSERT_EQ(0ull, qpos);
        }
    }
}

TEST(ExtendingPathsMatching, TypicalPath_NotExtendedWhenNonUniqMatch)
{
    {
        Graph graph = makeDeletionGraph("AAACC", "TTGGG", "TTAAA");
        const Path path(&graph, 4, { 0 }, 4);
        const string query = "CTT";
        size_t qpos = 0;
        const Path extended_path = extendPathMatching(path, query, qpos);
        const Path expected_path(&graph, 4, { 0 }, 5);
        ASSERT_EQ(expected_path, extended_path);
        ASSERT_EQ(0ull, qpos);
    }

    {
        Graph graph = makeDeletionGraph("AAACC", "ATGCC", "TTAAA");
        const Path path(&graph, 0, { 2 }, 0);
        const string query = "CCT";
        size_t qpos = 2;
        const Path extended_path = extendPathMatching(path, query, qpos);
        const Path expected_path(&graph, 0, { 2 }, 1);
        ASSERT_EQ(expected_path, extended_path);
        ASSERT_EQ(2ull, qpos);
    }

    {
        Graph graph = makeSwapGraph("AAAG", "AGCC", "A", "GTTT");
        const Path path(&graph, 0, { 0 }, 2);
        const string query = "AAAGAG";
        size_t qpos = 0;
        const Path extended_path = extendPathEndMatching(path, query, qpos);
        const Path expected_path(&graph, 0, { 0 }, 4);
        ASSERT_EQ(expected_path, extended_path);
        ASSERT_EQ(0ull, qpos); // extending to right doesn't move qpos
    }
}

TEST(SplittingSequenceByPath, SequenceOfDifferentLength_ExceptionRaised)
{
    Graph graph = makeDeletionGraph("AAAACC", "TTTGG", "ATTT");
    Path path(&graph, 3, { 0, 1 }, 2);
    const string sequence = "AA";
    EXPECT_ANY_THROW(splitSequenceByPath(path, sequence));
}

TEST(SplittingSequenceByPath, SingleNodePath_SequenceSplit)
{
    Graph graph = makeDeletionGraph("AAAACC", "TTTGG", "ATTT");
    Path path(&graph, 1, { 1 }, 4);
    const string sequence = "AAT";
    const vector<string> expected_pieces = { sequence };
    EXPECT_EQ(expected_pieces, splitSequenceByPath(path, sequence));
}

TEST(SplittingSequenceByPath, MultiNodePath_SequenceSplit)
{
    Graph graph = makeDeletionGraph("AAAACC", "TTTGG", "ATTT");
    {
        Path path(&graph, 1, { 0, 1 }, 4);
        const string sequence = "AAAAAGGGG";
        const vector<string> expected_pieces = { "AAAAA", "GGGG" };
        EXPECT_EQ(expected_pieces, splitSequenceByPath(path, sequence));
    }

    {
        Path path(&graph, 3, { 0, 2 }, 2);
        const string sequence = "AAACC";
        const vector<string> expected_pieces = { "AAA", "CC" };
        EXPECT_EQ(expected_pieces, splitSequenceByPath(path, sequence));
    }

    {
        Path path(&graph, 3, { 0, 1, 2 }, 2);
        const string sequence = "AAAGGGGGCC";
        const vector<string> expected_pieces = { "AAA", "GGGGG", "CC" };
        EXPECT_EQ(expected_pieces, splitSequenceByPath(path, sequence));
    }
}

TEST(GraphPathOperations, GraphPathsOverlapDetected)
{
    Graph swap = makeSwapGraph("AAAA", "TTTT", "CCCC", "GGGG");
    {
        const Path p1(&swap, 0, { 0, 1 }, 3);
        const Path p2(&swap, 0, { 1, 3 }, 3);

        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p2, p1));

        const Path expected_merge(&swap, 0, { 0, 1, 3 }, 3);
        ASSERT_EQ(mergePaths(p1, p2), expected_merge);
        ASSERT_EQ(mergePaths(p2, p1), expected_merge);
    }

    {
        const Path p1(&swap, 2, { 0, 1, 3 }, 2);
        const Path p2(&swap, 0, { 1, 3 }, 3);

        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p2, p1));

        const Path expected_merge(&swap, 2, { 0, 1, 3 }, 3);
        ASSERT_EQ(mergePaths(p1, p2), expected_merge);
        ASSERT_EQ(mergePaths(p2, p1), expected_merge);
    }
    {
        const Path p1(&swap, 2, { 0, 2 }, 1);
        const Path p2(&swap, 1, { 2 }, 3);

        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_TRUE(checkPathPrefixSuffixOverlap(p2, p1));

        const Path expected_merge(&swap, 2, { 0, 2 }, 3);
        ASSERT_EQ(mergePaths(p1, p2), expected_merge);
        ASSERT_EQ(mergePaths(p2, p1), expected_merge);
    }
}

TEST(GraphPathOperations, GraphPathsAdjacencyDetected)
{
    Graph graph = makeDoubleSwapGraph("AAAA", "TTTT", "CCCC", "GGGG", "TTTT", "CCCC", "AAAA");

    {
        // p1 ends just before p2 begins
        const Path p1(&graph, 0, { 0, 1 }, 1);
        const Path p2(&graph, 2, { 1, 3 }, 3);

        ASSERT_TRUE(checkIfPathsAdjacent(p1, p2));
        ASSERT_TRUE(checkIfPathsAdjacent(p2, p1));

        const Path expected_merge(&graph, 0, { 0, 1, 3 }, 3);
        ASSERT_EQ(mergePaths(p1, p2), expected_merge);
        ASSERT_EQ(mergePaths(p2, p1), expected_merge);
    }

    {
        // p1 ends too far before p2 begins
        const Path p1(&graph, 0, { 0, 1 }, 0);
        const Path p2(&graph, 2, { 1, 3 }, 3);

        ASSERT_FALSE(checkIfPathsAdjacent(p1, p2));
        ASSERT_FALSE(checkIfPathsAdjacent(p2, p1));
    }
    {
        // p1 ends just before p2 begins
        const Path p1(&graph, 0, { 0, 1 }, 3);
        const Path p2(&graph, 0, { 3 }, 3);

        ASSERT_TRUE(checkIfPathsAdjacent(p1, p2));
        ASSERT_TRUE(checkIfPathsAdjacent(p2, p1));
        const Path expected_merge(&graph, 0, { 0, 1, 3 }, 3);
        ASSERT_EQ(mergePaths(p1, p2), expected_merge);
        ASSERT_EQ(mergePaths(p2, p1), expected_merge);
    }

    {
        // p1 ends too far before p2 begins
        const Path p1(&graph, 0, { 0, 1 }, 2);
        const Path p2(&graph, 0, { 3 }, 3);

        ASSERT_FALSE(checkIfPathsAdjacent(p1, p2));
        ASSERT_FALSE(checkIfPathsAdjacent(p2, p1));
    }

    {
        // p1 ends too far before p2 begins
        const Path p1(&graph, 0, { 0, 1 }, 2);
        const Path p2(&graph, 0, { 4 }, 3);

        ASSERT_FALSE(checkIfPathsAdjacent(p1, p2));
        ASSERT_FALSE(checkIfPathsAdjacent(p2, p1));
    }
}

TEST(GraphPathOperations, GraphPathsNoOverlapDetected)
{
    Graph swap = makeSwapGraph("AAAA", "TTTT", "CCCC", "GGGG");
    {
        // p1 ends before p2 begins
        const Path p1(&swap, 0, { 0, 1 }, 1);
        const Path p2(&swap, 2, { 1, 3 }, 3);

        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p2, p1));
    }

    {
        // no shared nodes
        const Path p1(&swap, 0, { 0 }, 3);
        const Path p2(&swap, 2, { 1, 3 }, 3);

        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p2, p1));
    }

    {
        // incompatible
        const Path p1(&swap, 0, { 0, 1, 3 }, 3);
        const Path p2(&swap, 2, { 0, 2, 3 }, 3);

        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p2, p1));
    }
    {
        // incompatible 2
        const Path p1(&swap, 0, { 0, 1 }, 3);
        const Path p2(&swap, 2, { 2, 3 }, 3);

        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p1, p2));
        ASSERT_FALSE(checkPathPrefixSuffixOverlap(p2, p1));
    }
}

TEST(GraphPathOperations, PathsMergedExhaustively)
{
    Graph swap = makeDoubleSwapGraph("AAAA", "TTTT", "CCCC", "GGGG", "TTTT", "CCCC", "AAAA");

    {
        const Path p0(&swap, 0, { 1, 3 }, 3);
        const Path p1(&swap, 0, { 2, 3 }, 3);
        const Path p2(&swap, 0, { 3, 4 }, 3);
        const Path p3(&swap, 0, { 3, 5 }, 3);

        const list<Path> expected_merged_path{
            Path(&swap, 0, { 1, 3, 4 }, 3),
            Path(&swap, 0, { 2, 3, 5 }, 3),
            Path(&swap, 0, { 2, 3, 4 }, 3),
            Path(&swap, 0, { 1, 3, 5 }, 3),
        };
        list<Path> merged_path{ p0, p1, p2, p3 };
        graphtools::exhaustiveMerge(merged_path);

        ASSERT_EQ(merged_path, expected_merged_path);
    }
}

TEST(GraphPathOperations, IntersectPaths_NoIntersection)
{
    Graph swap = makeDoubleSwapGraph("AAAA", "TTTT", "CCCC", "GGGG", "TTTT", "CCCC", "AAAA");

    // no shared nodes
    {
        const Path p0(&swap, 0, { 1 }, 3);
        const Path p1(&swap, 0, { 2 }, 3);

        const auto intersection = graphtools::intersectPaths(p0, p1);
        ASSERT_TRUE(intersection.empty());

        const auto intersection_r = graphtools::intersectPaths(p1, p0);
        ASSERT_TRUE(intersection_r.empty());
    }

    // one shared node, but no shared sequence
    {
        const Path p0(&swap, 0, { 1, 3 }, 1);
        const Path p1(&swap, 2, { 3, 4 }, 3);

        const auto intersection = graphtools::intersectPaths(p0, p1);
        ASSERT_TRUE(intersection.empty());

        const auto intersection_r = graphtools::intersectPaths(p1, p0);
        ASSERT_TRUE(intersection_r.empty());
    }
}

TEST(GraphPathOperations, IntersectPaths_SimpleIntersection)
{
    Graph swap = makeDoubleSwapGraph("AAAA", "TTTT", "CCCC", "GGGG", "TTTT", "CCCC", "AAAA");

    // full node is shared
    {
        const Path p0(&swap, 0, { 1, 3, 5 }, 4);
        const Path p1(&swap, 0, { 2, 3, 4 }, 4);

        std::list<Path> expected{
            Path(&swap, 0, { 3 }, 4),

        };

        const auto intersection = graphtools::intersectPaths(p0, p1);
        const auto intersection_r = graphtools::intersectPaths(p1, p0);
        ASSERT_EQ(expected, intersection);
        ASSERT_EQ(expected, intersection_r);
    }

    // partial node is shared
    {
        const Path p0(&swap, 0, { 1, 3 }, 2);
        const Path p1(&swap, 1, { 3, 4 }, 3);

        std::list<Path> expected{
            Path(&swap, 1, { 3 }, 2),

        };

        const auto intersection = graphtools::intersectPaths(p0, p1);
        const auto intersection_r = graphtools::intersectPaths(p1, p0);
        ASSERT_EQ(expected, intersection);
        ASSERT_EQ(expected, intersection_r);
    }
}

TEST(GraphPathOperations, IntersectPaths_ComplexIntersection)
{
    Graph swap = makeDoubleSwapGraph("AAAA", "TTTT", "CCCC", "GGGG", "TTTT", "CCCC", "AAAA");

    // multiple full nodes shared, two resulting paths
    {
        const Path p0(&swap, 0, { 1, 3, 5, 6 }, 4);
        const Path p1(&swap, 0, { 2, 3, 4, 6 }, 4);

        std::list<Path> expected{
            Path(&swap, 0, { 3 }, 4),
            Path(&swap, 0, { 6 }, 4),
        };

        const auto intersection = graphtools::intersectPaths(p0, p1);
        const auto intersection_r = graphtools::intersectPaths(p1, p0);
        ASSERT_EQ(expected, intersection);
        ASSERT_EQ(expected, intersection_r);
    }

    // complex subpath match
    {
        const Path p0(&swap, 0, { 1, 3, 4 }, 2);
        const Path p1(&swap, 0, { 2, 3, 4, 6 }, 3);

        std::list<Path> expected{
            Path(&swap, 0, { 3, 4 }, 2),
        };

        const auto intersection = graphtools::intersectPaths(p0, p1);
        const auto intersection_r = graphtools::intersectPaths(p1, p0);
        ASSERT_EQ(expected, intersection);
        ASSERT_EQ(expected, intersection_r);
    }

    // complex subpath match
    {
        const Path p0(&swap, 0, { 1, 3, 4 }, 2);
        const Path p1(&swap, 2, { 3, 4, 6 }, 3);

        std::list<Path> expected{
            Path(&swap, 2, { 3, 4 }, 2),
        };

        const auto intersection = graphtools::intersectPaths(p0, p1);
        const auto intersection_r = graphtools::intersectPaths(p1, p0);
        ASSERT_EQ(expected, intersection);
        ASSERT_EQ(expected, intersection_r);
    }
}

TEST(GeneratingSubpathForEachNode, TypicalPaths_Split)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");

    {
        const Path path(&graph, 0, { 0 }, 1);
        const list<Path> expected_subpaths = { path };
        EXPECT_EQ(expected_subpaths, generateSubpathForEachNode(path));
    }

    {
        const Path path(&graph, 3, { 0, 1, 2 }, 0);
        const list<Path> expected_subpaths
            = { Path(&graph, 3, { 0 }, 3), Path(&graph, 0, { 1 }, 2), Path(&graph, 0, { 2 }, 0) };
        EXPECT_EQ(expected_subpaths, generateSubpathForEachNode(path));
    }

    {
        const Path path(&graph, 1, { 0, 1, 1, 1, 2 }, 2);
        const list<Path> expected_subpaths
            = { Path(&graph, 1, { 0 }, 3), Path(&graph, 0, { 1 }, 2), Path(&graph, 0, { 1 }, 2),
                Path(&graph, 0, { 1 }, 2), Path(&graph, 0, { 2 }, 2) };
        EXPECT_EQ(expected_subpaths, generateSubpathForEachNode(path));
    }
}

TEST(CheckingIfPathsAreBookended, AdjacentPathsWithEndsOnSameNode_CheckPassed)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    const Path first_path(&graph, 0, { 0, 1 }, 1);
    const Path second_path(&graph, 1, { 1, 2 }, 1);

    ASSERT_TRUE(checkIfBookended(first_path, second_path));
}

TEST(CheckingIfPathsAreBookended, AdjacentPathsThatEndOnDifferentNodes_CheckPassed)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    const Path first_path(&graph, 0, { 0 }, 3);
    const Path second_path(&graph, 0, { 1, 2 }, 1);

    ASSERT_TRUE(checkIfBookended(first_path, second_path));
}

TEST(CheckingIfPathsAreBookended, NonadjacentPathsWithEndsOnSameNode_CheckFailed)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    const Path first_path(&graph, 0, { 0, 1 }, 0);
    const Path second_path(&graph, 1, { 1, 2 }, 1);

    ASSERT_FALSE(checkIfBookended(first_path, second_path));
}

TEST(CheckingIfPathsAreBookended, NonadjacentPathsThatEndOnNeighboringNodes_CheckFailed)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    const Path first_path(&graph, 0, { 0 }, 2);
    const Path second_path(&graph, 0, { 1, 2 }, 1);

    ASSERT_FALSE(checkIfBookended(first_path, second_path));
}

TEST(CheckingIfPathsAreBookended, NonadjacentPathsThatEndOnNonneighboringNodes_CheckFailed)
{
    Graph graph = makeSwapGraph("TTT", "AT", "CAT", "CCCCC");
    const Path first_path(&graph, 0, { 1 }, 2);
    const Path second_path(&graph, 0, { 2 }, 3);

    ASSERT_FALSE(checkIfBookended(first_path, second_path));
}

TEST(MergingBookendedPaths, PathsThatAreNotBookended_ExceptionThrown)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    const Path first_path(&graph, 0, { 0 }, 2);
    const Path second_path(&graph, 0, { 1, 2 }, 1);

    EXPECT_ANY_THROW(concatenatePaths(first_path, second_path));
}

TEST(MergingBookendedPaths, AdjacentPathsWithEndsOnSameNode_Merged)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    const Path first_path(&graph, 0, { 0, 1 }, 1);
    const Path second_path(&graph, 1, { 1, 2 }, 1);

    Path merged_path = concatenatePaths(first_path, second_path);

    Path expected_path(&graph, 0, { 0, 1, 2 }, 1);
    EXPECT_EQ(expected_path, merged_path);
}

TEST(MergingBookendedPaths, AdjacentPathsThatEndOnDifferentNodes_Merged)
{
    Graph graph = makeStrGraph("TTT", "AT", "CCCCC");
    const Path first_path(&graph, 0, { 0 }, 3);
    const Path second_path(&graph, 0, { 1, 2 }, 1);

    Path merged_path = concatenatePaths(first_path, second_path);

    Path expected_path(&graph, 0, { 0, 1, 2 }, 1);
    EXPECT_EQ(expected_path, merged_path);
}
