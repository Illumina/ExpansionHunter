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

#include "graphalign/GraphAlignmentOperations.hh"

#include <string>

#include "gtest/gtest.h"

#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/Path.hh"

using std::list;
using std::string;

using namespace graphtools;

TEST(SplitingNodeCigarEncoding, TypicalCigarEncoding_CigarAndNodeIdExtracted)
{
    string cigar;
    NodeId node_id;
    splitNodeCigar("1[4M5S]", cigar, node_id);
    EXPECT_EQ(1ul, node_id);
    EXPECT_EQ("4M5S", cigar);
}

TEST(DecodingGraphAlignment, SinglenodeGraphAlignment_Decoded)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    GraphAlignment alignment = decodeGraphAlignment(1, "1[2M]", &graph);

    GraphAlignment expected_alignment(Path(&graph, 1, { 1 }, 3), { Alignment(1, "2M") });
    EXPECT_EQ(expected_alignment, alignment);
}

TEST(DecodingGraphAlignment, MultinodeGraphAlignment_Decoded)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    GraphAlignment alignment = decodeGraphAlignment(0, "0[4M]1[2M3S]", &graph);

    GraphAlignment expected_alignment(Path(&graph, 0, { 0, 1 }, 2), { Alignment(0, "4M"), Alignment(0, "2M3S") });
    EXPECT_EQ(expected_alignment, alignment);
}

TEST(CheckingConsistencyOfAlignments, ConsistentAlignment_CheckPassed)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    GraphAlignment alignment = decodeGraphAlignment(0, "0[4M]1[2M3S]", &graph);

    const string query = "AAAATTCCC";
    ASSERT_TRUE(checkConsistency(alignment, query));
}

TEST(CheckingConsistencyOfAlignments, QueryIsShorterThanAlignment_CheckFailed)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    GraphAlignment alignment = decodeGraphAlignment(0, "0[4M]1[2M3S]", &graph);

    const string query = "AAAATT";
    ASSERT_FALSE(checkConsistency(alignment, query));
}

TEST(CheckingConsistencyOfAlignments, GraphAlignmentWithInconsistentLinearAlignment_CheckFailed)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
    GraphAlignment alignment = decodeGraphAlignment(0, "0[4M]1[2M3S]", &graph);

    const string query = "AAAAGGCCC";
    ASSERT_FALSE(checkConsistency(alignment, query));
}

TEST(CheckingConsistencyOfAlignments, UnderlyingPathCanBeShortened_CheckFailed)
{
    Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");

    GraphAlignment alignment(Path(&graph, 4, { 0, 1 }, 2), { Alignment(4, "3S"), Alignment(0, "2M") });
    ASSERT_FALSE(checkConsistency(alignment, "GGGTT"));
}

TEST(GettingQuerySequencesForEachNode, TypicalAlignment_SequencePairs)
{
    //  AAA:CGG:CGG:TG
    //    |:   :|||:||
    // TTTA:C--:CGG:TGGTTT

    Graph graph = makeStrGraph("AAA", "CGG", "TG");
    GraphAlignment alignment = decodeGraphAlignment(2, "0[3S1M]1[1M2D]1[3M]2[2M4S]", &graph);

    const string query = "TTTACCGGTGGTTT";
    const list<string> query_pieces = getQuerySequencesForEachNode(alignment, query);

    const list<string> expected_pieces = { "TTTA", "C", "CGG", "TGGTTT" };
    ASSERT_EQ(expected_pieces, query_pieces);
}

TEST(PrettyPrintingAlignments, TypicalAlignment_PrettyPrinted)
{
    //  AAA:CGG:CGG:TG
    //    |:   :|||:||
    // TTTA:C--:CGG:TGGTTT

    Graph graph = makeStrGraph("AAA", "CGG", "TG");
    GraphAlignment alignment = decodeGraphAlignment(2, "0[3S1M]1[1M2D]1[3M]2[2M4S]", &graph);

    const string query = "TTTACCGGTGGTTT";
    const string encoding = prettyPrint(alignment, query);

    const string expected_encoding = "---A:CGG:CGG:TG----\n"
                                     "   |:|  :|||:||    \n"
                                     "TTTA:C--:CGG:TGGTTT";
    ASSERT_EQ(expected_encoding, encoding);
}

TEST(ProjectingLinearAlignmentOntoPath, TypicalLinearAlignment_GraphAlignment)
{
    //  CATAC
    //   || |
    // GGAT-CGAA
    //  00122
    const string query = "GGATCGAA";
    const string reference = "CAWAC";
    Alignment linear_alignment(1, "2S2M1D1M3S");

    Graph graph = makeStrGraph("AACA", "W", "ACTTT");
    Path path(&graph, 2, { 0, 1, 2 }, 1);

    GraphAlignment graph_alignment = projectAlignmentOntoGraph(linear_alignment, path);

    GraphAlignment expected_graph_alignment = decodeGraphAlignment(3, "0[2S1M]1[1M]2[1D1M3S]", &graph);

    ASSERT_EQ(expected_graph_alignment, graph_alignment);
}
