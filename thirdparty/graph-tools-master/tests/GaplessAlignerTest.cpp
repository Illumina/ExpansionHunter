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
