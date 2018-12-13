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
    GappedGraphAligner aligner(&graph, kmer_len, padding_len, seed_affix_trim_len, "path-aligner");

    // -> CGCGCGTA
    //    | ||||||
    // -> C-CGCGTA
    //    11111122

    Path seed_path(&graph, 3, { 0 }, 3);
    const size_t extension_len = 12;
    list<PathAndAlignment> extensions = aligner.extendAlignmentSuffix(seed_path, "CCGCGTA", extension_len);

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
    GappedGraphAligner aligner(&graph, kmer_len, padding_len, seed_affix_trim_len, "path-aligner");

    Path seed_path(&graph, 3, { 0 }, 3);
    list<PathAndAlignment> extensions = aligner.extendAlignmentSuffix(seed_path, "CCC", 3);

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

class AlignerTests : public ::testing::TestWithParam<std::string>
{
};

// TODO: Throw an error if there are no valid extensions?
TEST_P(AlignerTests, ExtendingAlignmentPrefix_TypicalSequences_AlignmentExtended)
{
    Graph graph = makeStrGraph("ATATTA", "CG", "TATTT");

    const size_t kmer_len = 3;
    const size_t padding_len = 5;
    const size_t seed_affix_trim_len = 0;
    GappedGraphAligner aligner(&graph, kmer_len, padding_len, seed_affix_trim_len, GetParam());

    //  ATTAC-GCGC <-
    //  || || |||
    //  ATAACAGCGG <-
    //  00001 1111

    Path seed_path(&graph, 1, { 1 }, 1);
    const size_t extension_len = 10;
    list<PathAndAlignment> extensions = aligner.extendAlignmentPrefix(seed_path, "ATAACAGCGG", extension_len);

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
    GappedGraphAligner aligner(
        &graph, kmer_len, padding_len, seed_affix_trim_len, GetParam(), LinearAlignmentParameters(5, -4, -8, 0));

    // TTA-CG-CG-TAT
    // ||  || |  |||
    // TT--CG-C--TAT

    list<GraphAlignment> alignments = aligner.align("TTCGCTAT");

    list<GraphAlignment> expected_alignments = { decodeGraphAlignment(3, "0[2M1D]1[2M]1[1M1D]2[3M]", &graph) };
    // with default m=5,mm=-4,go=-8,ge=-2 2M1D=10-8-2=0, 1M1X=5-4=1, so, test needs an update:
    // list<GraphAlignment> expected_alignments = { decodeGraphAlignment(4, "0[1M1X]1[2M]1[1M1D]2[3M]", &graph) };
    EXPECT_EQ(expected_alignments, alignments);
}

TEST_P(AlignerTests, PerformingGappedAlignment_MultimappingQuery_BestAlignmentsComputed)
{
    Graph graph = makeStrGraph("AAG", "CGG", "CTT");
    GappedGraphAligner aligner(&graph, 3, 0, 0, GetParam());

    // G-CG-C
    // 0-11-1
    // 1-11-1

    list<GraphAlignment> alignments = aligner.align("GCGGC");

    list<GraphAlignment> expected_alignments
        = { decodeGraphAlignment(2, "0[1M]1[3M]1[1M]", &graph), decodeGraphAlignment(2, "0[1M]1[3M]2[1M]", &graph),
            decodeGraphAlignment(2, "1[1M]1[3M]1[1M]", &graph), decodeGraphAlignment(2, "1[1M]1[3M]2[1M]", &graph) };
    EXPECT_EQ(expected_alignments, alignments);
}

TEST_P(AlignerTests, PerformingGappedAlignment_KmerExtensionInBothDirectionsNotNeeded_BestAlignmentsComputed)
{
    Graph graph = makeStrGraph("AAG", "CGG", "CTT");
    GappedGraphAligner aligner(&graph, 3, 0, 0, GetParam());

    {
        list<GraphAlignment> alignments = aligner.align("CGGCT");
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[3M]2[2M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("AATCGG");
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "0[2M1X]1[3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("CTT");
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "2[3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }
}

TEST_P(AlignerTests, PerformingGappedAlignment_KmerExtensionIsUnalignable_BestAlignmentsComputed)
{
    Graph graph = makeStrGraph("AAG", "CGG", "CTT");
    GappedGraphAligner aligner(&graph, 3, 0, 0, GetParam());

    {
        list<GraphAlignment> alignments = aligner.align("CGGAA");
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[3M2S]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("TTCGG");
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[2S3M]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("TCGGA");
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[1S3M1S]", &graph) };
        EXPECT_EQ(expected_alignments, alignments);
    }
}

TEST_P(AlignerTests, PerformingGappedAlignment_PolyalanineRepeat_ReadAligned)
{
    Graph graph = makeStrGraph("AAG", "GCN", "ATT");
    GappedGraphAligner aligner(&graph, 4, 0, 0, GetParam());

    list<GraphAlignment> alignments = aligner.align("AGGCCGTGGCAATT");
    list<GraphAlignment> expected_alignments = { decodeGraphAlignment(1, "0[2M]1[3M]1[1M1X1M]1[3M]2[3M]", &graph) };

    EXPECT_EQ(expected_alignments, alignments);
}

TEST_P(AlignerTests, PerformingGappedAlignment_ReadWithLowqualityBases_ReadAligned)
{
    Graph graph = makeStrGraph("AAG", "CGG", "CTT");
    GappedGraphAligner aligner(&graph, 4, 0, 0, GetParam());

    list<GraphAlignment> alignments = aligner.align("aagcggctt");
    list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "0[3M]1[3M]2[3M]", &graph) };

    EXPECT_EQ(expected_alignments, alignments);
}

TEST_P(AlignerTests, PerformingGappedAlignment_IncorrectSeedKmer_ReadAligned)
{
    Graph graph = makeStrGraph("AAAA", "CCG", "TTTT");
    const int32_t seed_affix_trim_len = 2;
    GappedGraphAligner aligner(&graph, 4, 0, seed_affix_trim_len, GetParam());

    {
        list<GraphAlignment> alignments = aligner.align("CCACCGTTTT");
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[2M1X]1[3M]2[4M]", &graph) };

        EXPECT_EQ(expected_alignments, alignments);
    }

    {
        list<GraphAlignment> alignments = aligner.align("CCGTCG");
        list<GraphAlignment> expected_alignments = { decodeGraphAlignment(0, "1[3M]1[1X2M]", &graph) };

        EXPECT_EQ(expected_alignments, alignments);
    }
}

INSTANTIATE_TEST_CASE_P(
    AlignerTestsInst, AlignerTests, ::testing::Values(std::string("path-aligner"), std::string("dag-aligner")), );
