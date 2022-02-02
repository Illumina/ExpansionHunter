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

#include "graphalign/GappedAligner.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignment.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"
#include "graphcore/Path.hh"

using std::list;
using std::make_pair;
using std::string;

using namespace graphtools;

TEST(DISABLED_ExtendingAlignmentSuffix, UniquelyMappingQuery_AlignmentExtended)
{
    Graph graph = makeStrGraph("ATA", "CG", "TATTTTTTTTT");

    const size_t kmer_len = 3;
    const size_t padding_len = 5;
    const int32_t seed_affix_trim_len = 0;
    GappedGraphAligner aligner(&graph, kmer_len, padding_len, seed_affix_trim_len);
    AlignerSelector alignerSelector(AlignerType::PATH_ALIGNER);

    // -> CGCGCGTA
    //    | ||||||
    // -> C-CGCGTA
    //    11111122

    Path seed_path(&graph, 3, { 0 }, 3);
    const size_t extension_len = 12;
    list<PathAndAlignment> extensions
        = aligner.extendAlignmentSuffix(seed_path, "CCGCGTA", extension_len, alignerSelector);

    Alignment expected_alignment(0, "1M1D6M");
    Path expected_path(&graph, 3, { 0, 1, 1, 1, 2 }, 2);
    list<PathAndAlignment> expected_extensions = { make_pair(expected_path, expected_alignment) };

    EXPECT_EQ(expected_extensions, extensions);
}

TEST(DISABLED_ExtendingAlignmentSuffix, MultiMappingQuery_AlignmentExtended)
{
    Graph graph = makeStrGraph("AAA", "C", "CCA");

    const size_t kmer_len = 3;
    const size_t padding_len = 0;
    const int32_t seed_affix_trim_len = 0;
    GappedGraphAligner aligner(&graph, kmer_len, padding_len, seed_affix_trim_len);
    AlignerSelector alignerSelector(AlignerType::PATH_ALIGNER);

    Path seed_path(&graph, 3, { 0 }, 3);
    list<PathAndAlignment> extensions = aligner.extendAlignmentSuffix(seed_path, "CCC", 3, alignerSelector);

    Alignment expected_alignment_a = Alignment(0, "3M");
    Path expected_path_a(&graph, 3, { 0, 1, 1, 1 }, 1);

    Alignment expected_alignment_b = Alignment(0, "3M");
    Path expected_path_b(&graph, 3, { 0, 1, 1, 2 }, 1);

    Alignment expected_alignment_c = Alignment(0, "3M");
    Path expected_path_c(&graph, 3, { 0, 1, 2 }, 2);

    list<PathAndAlignment> expected_extensions
        = { make_pair(expected_path_a, expected_alignment_a), make_pair(expected_path_b, expected_alignment_b),
            make_pair(expected_path_c, expected_alignment_c) };

    EXPECT_EQ(expected_extensions, extensions);
}

class AlignerTests : public ::testing::TestWithParam<graphtools::AlignerType>
{
};

// TODO: Throw an error if there are no valid extensions?
TEST_P(AlignerTests, ExtendingAlignmentPrefix_TypicalSequences_AlignmentExtended)
{
    Graph graph = makeStrGraph("ATATTA", "CG", "TATTT");

    const size_t kmer_len = 3;
    const size_t padding_len = 5;
    const size_t seed_affix_trim_len = 0;
    GappedGraphAligner aligner(&graph, kmer_len, padding_len, seed_affix_trim_len);
    AlignerSelector alignerSelector(GetParam());

    //  ATTAC-GCGC <-
    //  || || |||
    //  ATAACAGCGG <-
    //  00001 1111

    Path seed_path(&graph, 1, { 1 }, 1);
    const size_t extension_len = 10;
    list<PathAndAlignment> extensions
        = aligner.extendAlignmentPrefix(seed_path, "ATAACAGCGG", extension_len, alignerSelector);

    Alignment expected_alignment(0, "2M1X2M1I3M1X");
    Path expected_path(&graph, 2, { 0, 1, 1, 1 }, 1);
    list<PathAndAlignment> expected_extensions = { make_pair(expected_path, expected_alignment) };

    EXPECT_EQ(expected_extensions, extensions);
}

TEST_P(AlignerTests, PerformingGappedAlignment_UniquelyMappingQuery_AlignmentPerformed)
{
    Graph graph = makeStrGraph("ATATTA", "CG", "TATTT");

    const size_t kmer_len = 3;
    const size_t padding_len = 2;
    const int32_t seed_affix_trim_len = 0;
    GappedGraphAligner aligner(&graph, kmer_len, padding_len, seed_affix_trim_len);
    AlignerSelector alignerSelector(GetParam(), LinearAlignmentParameters(5, -4, -8, 0));

    {
        // TTA-CG-CG-TAT
        // ||  || |  |||
        // TT--CG-C--TAT

        list<GraphAlignment> alignments = aligner.align("TTCGCTAT", alignerSelector);

        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(3, "0[2M1D]1[2M]1[1M1D]2[3M]", &graph) };
        // with default m=5,mm=-4,go=-8,ge=-2 2M1D=10-8-2=0, 1M1X=5-4=1, so, test needs an update:
        // list<GraphAlignment> expected_alignments = { decodeGraphAlignment(4, "0[1M1X]1[2M]1[1M1D]2[3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("ATTCGCTAT", alignerSelector);

        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(2, "0[3M1D]1[2M]1[1M1D]2[3M]", &graph) };
        // with default m=5,mm=-4,go=-8,ge=-2 2M1D=10-8-2=0, 1M1X=5-4=1, so, test needs an update:
        // list<GraphAlignment> expected_alignments = { decodeGraphAlignment(4, "0[1M1X]1[2M]1[1M1D]2[3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }
}

TEST_P(AlignerTests, PerformingGappedAlignment_MultimappingQuery_BestAlignmentsComputed)
{
    Graph graph = makeStrGraph("AAG", "CGG", "CTT");
    GappedGraphAligner aligner(&graph, 3, 0, 0);
    AlignerSelector alignerSelector(GetParam());

    // G-CG-C
    // 0-11-1
    // 1-11-1

    list<GraphAlignment> alignments = aligner.align("GCGGC", alignerSelector);

    list<GraphAlignment> expected_alignments
        = { decodeGraphAlignment(2, "0[1M]1[3M]1[1M]", &graph), decodeGraphAlignment(2, "0[1M]1[3M]2[1M]", &graph),
            decodeGraphAlignment(2, "1[1M]1[3M]1[1M]", &graph), decodeGraphAlignment(2, "1[1M]1[3M]2[1M]", &graph) };
    EXPECT_EQ(expected_alignments, alignments);
}

TEST_P(AlignerTests, PerformingGappedAlignment_KmerExtensionInBothDirectionsNotNeeded_BestAlignmentsComputed)
{
    Graph graph = makeStrGraph("AAG", "CGG", "CTT");
    GappedGraphAligner aligner(&graph, 3, 0, 0);
    AlignerSelector alignerSelector(GetParam());

    {
        list<GraphAlignment> alignments = aligner.align("CGGCT", alignerSelector);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[3M]2[2M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("AATCGG", alignerSelector);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "0[2M1X]1[3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("CTT", alignerSelector);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "2[3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }
}

TEST_P(AlignerTests, PerformingGappedAlignment_KmerExtensionIsUnalignable_BestAlignmentsComputed)
{
    Graph graph = makeStrGraph("AAG", "CGG", "CTT");
    GappedGraphAligner aligner(&graph, 3, 0, 0);
    AlignerSelector alignerSelector(GetParam());

    {
        list<GraphAlignment> alignments = aligner.align("CGGAA", alignerSelector);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[3M2S]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("TTCGG", alignerSelector);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[2S3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("TCGGA", alignerSelector);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[1S3M1S]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }
}

TEST_P(AlignerTests, PerformingGappedAlignment_PolyalanineRepeat_ReadAligned)
{
    Graph graph = makeStrGraph("AAG", "GCN", "ATT");
    GappedGraphAligner aligner(&graph, 4, 0, 0);
    AlignerSelector alignerSelector(GetParam());

    list<GraphAlignment> alignments = aligner.align("AGGCCGTGGCAATT", alignerSelector);
    list<GraphAlignment> expected_alignments = { decodeGraphAlignment(1, "0[2M]1[3M]1[1M1X1M]1[3M]2[3M]", &graph) };

    EXPECT_EQ(expected_alignments, alignments);
}

TEST_P(AlignerTests, PerformingGappedAlignment_ReadWithLowqualityBases_ReadAligned)
{
    Graph graph = makeStrGraph("AAG", "CGG", "CTT");
    GappedGraphAligner aligner(&graph, 4, 0, 0);
    AlignerSelector alignerSelector(GetParam());

    list<GraphAlignment> alignments = aligner.align("aagcggctt", alignerSelector);
    list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "0[3M]1[3M]2[3M]", &graph) };

    EXPECT_EQ(expected_alignments, alignments);
}

TEST_P(AlignerTests, PerformingGappedAlignment_IncorrectSeedKmer_ReadAligned)
{
    Graph graph = makeStrGraph("AAAA", "CCG", "TTTT");
    const int32_t seed_affix_trim_len = 2;
    GappedGraphAligner aligner(&graph, 4, 0, seed_affix_trim_len);
    AlignerSelector alignerSelector(GetParam());

    {
        list<GraphAlignment> alignments = aligner.align("CCACCGTTTT", alignerSelector);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[2M1X]1[3M]2[4M]", &graph) };

        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("CCGTCG", alignerSelector);
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[3M]1[1X2M]", &graph) };

        EXPECT_EQ(expected_alignments, alignments);
    }
}

TEST_P(AlignerTests, PerformingGappedAlignment_NoExceptionThrown)
{
    //     0       1          2          3  4     5       6     7
    //<left flank>(AT)*(GATATATATATATAT)*G(AT)*TTATATATG(AT)*<right flank>
    Graph graph(8);
    graph.setNodeSeq(0, "AGGATGACAGTAATATTATCTTACTATCTTACTATGTGTTACTTTATTAGTTTTTCCCTTATATGTTTGTTTTGGGATATATGACTTGGCTC");
    graph.setNodeSeq(1, "AT");
    graph.setNodeSeq(2, "GATATATATATATAT");
    graph.setNodeSeq(3, "G");
    graph.setNodeSeq(4, "AT");
    graph.setNodeSeq(5, "TTATATATG");
    graph.setNodeSeq(6, "AT");
    graph.setNodeSeq(7, "GATATATATTTATATTAAAAGGTGCTTTGTTCTTTGCAAAACAGTCTCCTATGTTATTTCCTCATTTTATTAAAATGTAACCTAAAACTGTT");

    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    graph.addEdge(2, 3);
    graph.addEdge(3, 4);
    graph.addEdge(4, 5);
    graph.addEdge(5, 6);
    graph.addEdge(6, 7);

    graph.addEdge(1, 1);
    graph.addEdge(0, 2);
    graph.addEdge(0, 3);

    graph.addEdge(2, 2);
    graph.addEdge(1, 3);

    graph.addEdge(4, 4);
    graph.addEdge(3, 5);

    graph.addEdge(6, 6);
    graph.addEdge(5, 7);

    const size_t kmer_len = 14;
    const size_t padding_len = 10;
    const int32_t seed_affix_trim_len = 14;
    GappedGraphAligner aligner(&graph, kmer_len, padding_len, seed_affix_trim_len);
    AlignerSelector alignerSelector(GetParam());

    const string query = "ctTTttgaTTTtttccctcacatgTTTTTtatatGataTtTctcTtCtCtcataTAtttaTAtAtAttAtATtTAtAtataTctttTAtATAT"
                         "AtaATaTaTaTATatCATATAtATaTATGATATATATATATATCATATATATATATG";

    ASSERT_NO_THROW(aligner.align(query, alignerSelector));
}

TEST_P(AlignerTests, PerformingGappedAlignment_FlankWithStrKmer_ReadAligned)
{
    Graph graph = makeStrGraph("AAAA", "CGG", "TTCGGCGGTT");
    const int32_t seed_affix_trim_len = 2;
    GappedGraphAligner aligner(&graph, 4, 0, seed_affix_trim_len);
    AlignerSelector alignerSelector(GetParam());

    list<GraphAlignment> alignments = aligner.align("CGGCGGCGGCGGCGG", alignerSelector);
    list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[3M]1[3M]1[3M]1[3M]1[3M]", &graph) };

    EXPECT_EQ(expected_alignments, alignments);
}

INSTANTIATE_TEST_SUITE_P(
    AlignerTestsInst, AlignerTests, ::testing::Values(AlignerType::PATH_ALIGNER, AlignerType::DAG_ALIGNER));
