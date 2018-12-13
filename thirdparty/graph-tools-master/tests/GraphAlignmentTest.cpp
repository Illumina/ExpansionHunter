//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
