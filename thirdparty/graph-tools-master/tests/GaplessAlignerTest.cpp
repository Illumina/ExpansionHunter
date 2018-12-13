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

#include "graphalign/GaplessAligner.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphutils/SequenceOperations.hh"

#include "gtest/gtest.h"

using std::list;
using std::string;

using namespace graphtools;

TEST(AligningTwoSequences, SequencesWithUnequalLength_ExceptionThrown)
{
    EXPECT_ANY_THROW(alignWithoutGaps(0, "AAA", "AAAA"));
}

TEST(AligningTwoSequences, EmptySequences_ExceptionThrown) { EXPECT_ANY_THROW(alignWithoutGaps(0, "", "")); }

TEST(AligningSequences, TypicalSequences_Aligned)
{
    const string reference = "NNNNATCGTTTG";
    const string query = "AGGTTTTG";
    const Alignment expected_alignment(4, "1M3X4M");
    ASSERT_EQ(expected_alignment, alignWithoutGaps(4, reference, query));
}

TEST(AligningSequences, SequenceWithDegenerateBases_Aligned)
{
    const string reference = "VVVVV";
    const string query = "AATTC";
    const Alignment expected_alignment(0, "2M2X1M");
    ASSERT_EQ(expected_alignment, alignWithoutGaps(0, reference, query));
}

TEST(AligningSequenceToPath, SingleNodePath_Aligned)
{
    Graph graph = makeDeletionGraph("AAAACC", "TTTGG", "ATTT");
    Path path(&graph, 1, { 1 }, 5);
    const string query = "ATGC";

    GraphAlignment expected_graph_alignment = decodeGraphAlignment(1, "1[1X2M1X]", &graph);
    GraphAlignment graph_alignment = alignWithoutGaps(path, query);
    EXPECT_EQ(expected_graph_alignment, graph_alignment);
}

TEST(AligningSequenceToPath, MultiNodePath_Aligned)
{
    Graph graph = makeDeletionGraph("AAAACC", "TTTGG", "ATTT");
    Path path(&graph, 2, { 0, 1, 2 }, 2);
    const string query = "TTCCTTAGGAT";

    GraphAlignment expected_graph_alignment = decodeGraphAlignment(2, "0[2X2M]1[2M1X2M]2[2M]", &graph);
    GraphAlignment graph_alignment = alignWithoutGaps(path, query);
    EXPECT_EQ(expected_graph_alignment, graph_alignment);
}

TEST(AligningSequenceToPath, TypicalStrPath_Aligned)
{
    Graph graph = makeStrGraph("AAAACC", "CCG", "ATTT");
    Path path(&graph, 2, { 0, 1, 1, 1, 2 }, 4);
    //                   FFFFRRRRRRRRRFFFF
    const string query = "AACCCCGCCGCCGATTT";

    GraphAlignment expected_graph_alignment = decodeGraphAlignment(2, "0[4M]1[3M]1[3M]1[3M]2[4M]", &graph);
    GraphAlignment graph_alignment = alignWithoutGaps(path, query);
    EXPECT_EQ(expected_graph_alignment, graph_alignment);
}

TEST(KmerExtraction, TypicalSequence_KmersExtracted)
{
    const string sequence = "AAatTT";
    const list<string> expected_4mers = { "AAAT", "AATT", "ATTT" };
    ASSERT_EQ(expected_4mers, extractKmersFromAllPositions(sequence, 4));

    const list<string> expected_7mers = {};
    ASSERT_EQ(expected_7mers, extractKmersFromAllPositions(sequence, 7));
}

TEST(AlignmentOfSequenceToShortPath, TypicalSequence_BestAlignmentObtained)
{
    Graph graph = makeDeletionGraph("AAACC", "TTGGG", "TTAAA");
    const Path path(&graph, 4, { 0 }, 4);
    const string query = "CCTTA";

    list<GraphAlignment> alignments = getBestAlignmentToShortPath(path, 1, query);

    list<GraphAlignment> expected_alignments = { decodeGraphAlignment(3, "0[2M]2[3M]", &graph) };
    ASSERT_EQ(expected_alignments, alignments);
}

TEST(AlignmentOfSequenceToGraph, TypicalSequence_BestAlignmentObtained)
{
    Graph graph = makeDeletionGraph("AAAACC", "TTTGG", "ATTT");

    const int32_t kmer_len = 3;
    GaplessAligner aligner(&graph, kmer_len);

    const string query = "TTCCTTAGGAT";
    list<GraphAlignment> alignments = aligner.align(query);

    list<GraphAlignment> expected_alignments = { decodeGraphAlignment(2, "0[2X2M]1[2M1X2M]2[2M]", &graph) };
    ASSERT_EQ(expected_alignments, alignments);
}

TEST(GraphAlignment, TypicalStrGraph_BestAlignmentObtained)
{
    Graph graph = makeStrGraph("AAAACG", "CCG", "ATTT");
    const int32_t kmer_len = 3;
    GaplessAligner aligner(&graph, kmer_len);

    {
        //                            FFFFRRRRRRRRRFFFF
        const string spanning_read = "AACGCCGCCGCCGATTT";
        list<GraphAlignment> alignments = aligner.align(spanning_read);

        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(2, "0[4M]1[3M]1[3M]1[3M]2[4M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        //                          RRRRRRRRRRR
        const string repeat_read = "CGCCGCCGCCG";
        list<GraphAlignment> alignments = aligner.align(repeat_read);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(4, "0[2M]1[3M]1[3M]1[3M]", &graph),
                                                     decodeGraphAlignment(1, "1[2M]1[3M]1[3M]1[3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        //                          RRRXRRRRXRRR
        const string repeat_read = "CCGACGCCTCCG";
        list<GraphAlignment> alignments = aligner.align(repeat_read);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[3M]1[1X2M]1[2M1X]1[3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }
}

TEST(GraphAlignment, PolyalanineGraph_BestAlignmentObtained)
{
    Graph graph = makeStrGraph("AACG", "GCN", "ATTT");
    const int32_t kmer_len = 3;
    GaplessAligner aligner(&graph, kmer_len);

    {
        //                            FFFFGCNGCNGCNFFFF
        const string spanning_read = "AACGGCAGCTGCGATTT";
        list<GraphAlignment> alignments = aligner.align(spanning_read);

        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "0[4M]1[3M]1[3M]1[3M]2[4M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        //                          CNGCNGCNGCN
        const string repeat_read = "CGGCAGCTGCG";
        list<GraphAlignment> alignments = aligner.align(repeat_read);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(2, "0[2M]1[3M]1[3M]1[3M]", &graph),
                                                     decodeGraphAlignment(1, "1[2M]1[3M]1[3M]1[3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }
}
