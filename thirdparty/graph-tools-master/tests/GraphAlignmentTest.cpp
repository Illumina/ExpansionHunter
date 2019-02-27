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

#include "graphalign/GraphAlignment.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/Path.hh"

using std::list;
using std::string;
using std::vector;

using namespace graphtools;

TEST(InitializingGraphAlignment, CompatiblePath_GraphAlignmentCreated)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");

    {
        Path path(&graph, 3, { 0, 1, 2 }, 3);
        vector<Alignment> alignments = { Alignment(3, "1M"), Alignment(0, "4M"), Alignment(0, "3M") };
        EXPECT_NO_THROW(GraphAlignment(path, alignments));
    }

    {
        Path path(&graph, 2, { 1 }, 3);
        vector<Alignment> alignments = { Alignment(2, "1M") };
        EXPECT_NO_THROW(GraphAlignment(path, alignments));
    }
}

TEST(InitializingGraphAlignment, IncompatiblePath_ExceptionThrown)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    Path path(&graph, 2, { 0, 1, 2 }, 3);

    {
        vector<Alignment> alignments = { Alignment(3, "1M"), Alignment(0, "4M"), Alignment(0, "3M") };
        EXPECT_ANY_THROW(GraphAlignment(path, alignments));
    }

    {
        vector<Alignment> alignments = { Alignment(2, "2M"), Alignment(0, "4M"), Alignment(0, "4M") };
        EXPECT_ANY_THROW(GraphAlignment(path, alignments));
    }

    {
        vector<Alignment> alignments = { Alignment(2, "2M"), Alignment(0, "3M"), Alignment(0, "3M") };
        EXPECT_ANY_THROW(GraphAlignment(path, alignments));
    }

    {
        vector<Alignment> alignments = { Alignment(2, "2M"), Alignment(1, "4M"), Alignment(0, "3M") };
        EXPECT_ANY_THROW(GraphAlignment(path, alignments));
    }
}

TEST(GettingNumMatchesInGraphAlignment, TypicalGraphAlignment_GotNumMatches)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    const string query = "AAAATTCCC";
    GraphAlignment graph_alignment = decodeGraphAlignment(0, "0[4M]1[2M3S]", &graph);
    EXPECT_EQ(6u, graph_alignment.numMatches());
}

TEST(GettingGraphAlignmentSpans, TypicalGraphAlignment_GotQueryAndReferenceSpans)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    const string query = "AAAATTCCC";
    GraphAlignment graph_alignment = decodeGraphAlignment(0, "0[4M]1[2M3S]", &graph);
    EXPECT_EQ(9u, graph_alignment.queryLength());
    EXPECT_EQ(6u, graph_alignment.referenceLength());
}

TEST(AccessingNodeAlignmentsByIndex, TypicalGraphAlignment_NodeAlignmentsAccessed)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGC", "TTTT");
    const string query = "AAAATTCCC";
    GraphAlignment graph_alignment = decodeGraphAlignment(0, "0[4M]1[2M3S]", &graph);
    EXPECT_EQ(Alignment(0, "4M"), graph_alignment[0]);
    EXPECT_EQ(Alignment(0, "2M3S"), graph_alignment[1]);
}

TEST(GettingIndexesOfNode, TypicalAlignment_IndexesObtained)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    const string read = "CCCCGCCGAT";
    GraphAlignment alignment = decodeGraphAlignment(4, "0[2M]1[3M]1[3M]2[2M]", &graph);
    const list<int32_t> left_flank_indexes = { 0 };
    const list<int32_t> repeat_unit_indexes = { 1, 2 };
    const list<int32_t> right_flank_indexes = { 3 };
    EXPECT_EQ(left_flank_indexes, alignment.getIndexesOfNode(0));
    EXPECT_EQ(repeat_unit_indexes, alignment.getIndexesOfNode(1));
    EXPECT_EQ(right_flank_indexes, alignment.getIndexesOfNode(2));
}

TEST(GettingIndexesOfNode, NodeNotInAlignment_EmptyListReturned)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    const string read = "ACCCCG";
    GraphAlignment alignment = decodeGraphAlignment(3, "0[3M]1[3M]", &graph);
    const list<int32_t> empty_list;
    EXPECT_EQ(empty_list, alignment.getIndexesOfNode(2));
    EXPECT_EQ(empty_list, alignment.getIndexesOfNode(4));
}

TEST(CheckingIfAlignmentOverlapsNode, TypicalAlignment_ChecksPerformed)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    const string read = "ACCCCG";
    GraphAlignment alignment = decodeGraphAlignment(3, "0[3M]1[3M]", &graph);
    EXPECT_TRUE(alignment.overlapsNode(0));
    EXPECT_TRUE(alignment.overlapsNode(1));
    EXPECT_FALSE(alignment.overlapsNode(2));
    EXPECT_FALSE(alignment.overlapsNode(3));
}

TEST(EncodingGraphAlignment, TypicalGraphAlignment_CigarStringObtained)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    const string read = "CCCCGCCGAT";
    const string cigar_string = "0[2M]1[3M]1[3M]2[2M]";
    GraphAlignment alignment = decodeGraphAlignment(4, cigar_string, &graph);

    ASSERT_EQ(cigar_string, alignment.generateCigar());
}

TEST(ComparingGraphAlignments, TypicalGraphAlignments_Compared)
{
    Graph graph = makeStrGraph("ATT", "CCG", "CTTT");

    GraphAlignment alignment_a = decodeGraphAlignment(1, "0[2M]1[3M]1[1M]", &graph);
    GraphAlignment alignment_b = decodeGraphAlignment(1, "0[2M]1[3M]2[1M]", &graph);

    EXPECT_TRUE(alignment_a < alignment_b);
    EXPECT_FALSE(alignment_b < alignment_a);
    EXPECT_FALSE(alignment_a == alignment_b);
}

TEST(ShrinkingGraphAlignmentStarts, TypicalGraphAlignment_Shrank)
{
    Graph graph = makeStrGraph("ATT", "CCG", "CTTT");

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[2M]1[3M]1[1M]", &graph);
        alignment.shrinkStart(1);
        GraphAlignment expectedAlignment = decodeGraphAlignment(2, "0[1S1M]1[3M]1[1M]", &graph);

        EXPECT_EQ(alignment, expectedAlignment);
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[2M]1[3M]1[1M]", &graph);
        alignment.shrinkStart(2);
        GraphAlignment expectedAlignment = decodeGraphAlignment(0, "1[2S3M]1[1M]", &graph);

        EXPECT_EQ(alignment, expectedAlignment);
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[2M]1[3M]1[1M]", &graph);
        alignment.shrinkStart(5);
        GraphAlignment expectedAlignment = decodeGraphAlignment(0, "1[5S1M]", &graph);

        EXPECT_EQ(alignment, expectedAlignment);
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[1S2M]1[3M]1[1M]", &graph);
        alignment.shrinkStart(3);
        GraphAlignment expectedAlignment = decodeGraphAlignment(1, "1[4S2M]1[1M]", &graph);

        EXPECT_EQ(alignment, expectedAlignment);
    }
}

TEST(ShrinkingGraphAlignmentStarts, ShrinkingByAlignmentLengthOrMore_ExceptionThrown)
{
    Graph graph = makeStrGraph("ATT", "CCG", "CTTT");

    GraphAlignment alignment = decodeGraphAlignment(1, "0[2M]1[3M]1[1M]", &graph);
    EXPECT_ANY_THROW(alignment.shrinkStart(alignment.referenceLength()));
    EXPECT_ANY_THROW(alignment.shrinkStart(alignment.referenceLength() + 1));
}

TEST(ShrinkingGraphAlignmentEnds, TypicalGraphAlignment_Shrank)
{
    Graph graph = makeStrGraph("ATT", "CCG", "CTTT");

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[2M]1[3M]1[1M]", &graph);
        alignment.shrinkEnd(1);
        GraphAlignment expectedAlignment = decodeGraphAlignment(1, "0[2M]1[3M1S]", &graph);

        EXPECT_EQ(alignment, expectedAlignment);
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[2M]1[3M]1[1M]", &graph);
        alignment.shrinkEnd(2);
        GraphAlignment expectedAlignment = decodeGraphAlignment(1, "0[2M]1[2M2S]", &graph);

        EXPECT_EQ(alignment, expectedAlignment);
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[2M]1[3M]1[1M]", &graph);
        alignment.shrinkEnd(5);
        GraphAlignment expectedAlignment = decodeGraphAlignment(1, "0[1M5S]", &graph);

        EXPECT_EQ(alignment, expectedAlignment);
    }

    {
        GraphAlignment alignment = decodeGraphAlignment(1, "0[1S2M]1[3M]1[1M3S]", &graph);
        alignment.shrinkEnd(4);
        GraphAlignment expectedAlignment = decodeGraphAlignment(1, "0[1S2M7S]", &graph);

        EXPECT_EQ(alignment, expectedAlignment);
    }
}

TEST(ShrinkingGraphAlignmentEnds, ShrinkingByAlignmentLengthOrMore_ExceptionThrown)
{
    Graph graph = makeStrGraph("ATT", "CCG", "CTTT");

    GraphAlignment alignment = decodeGraphAlignment(1, "0[2M]1[3M]1[1M]", &graph);
    EXPECT_ANY_THROW(alignment.shrinkEnd(alignment.referenceLength()));
    EXPECT_ANY_THROW(alignment.shrinkEnd(alignment.referenceLength() + 1));
}
