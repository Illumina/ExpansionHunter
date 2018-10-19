// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

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

#include "graphalign/KmerIndex.hh"

#include <list>

#include "gtest/gtest.h"

#include "graphcore/GraphBuilders.hh"
#include "graphcore/Path.hh"

using std::list;
using std::string;
using std::unordered_set;

using namespace graphtools;

TEST(KmerIndexInitialization, 1mers_IndexCreated)
{
    Graph graph = makeDeletionGraph("AC", "GG", "CAG");

    const int32_t kmer_size = 1;
    KmerIndex kmer_index(graph, kmer_size);

    const list<Path> a_paths = { Path(&graph, 0, { 0 }, 1), Path(&graph, 1, { 2 }, 2) };
    const list<Path> c_paths = { Path(&graph, 1, { 0 }, 2), Path(&graph, 0, { 2 }, 1) };
    const list<Path> g_paths = { Path(&graph, 0, { 1 }, 1), Path(&graph, 1, { 1 }, 2), Path(&graph, 2, { 2 }, 3) };

    const StringToPathsMap kmer_to_paths_maps = { { "A", a_paths }, { "C", c_paths }, { "G", g_paths } };

    KmerIndex expected_kmer_index(kmer_to_paths_maps);
    ASSERT_EQ(expected_kmer_index, kmer_index);
}

TEST(KmerIndexInitialization, 2mers_IndexCreated)
{
    Graph graph = makeDeletionGraph("AK", "GG", "CAG");
    const int32_t kmer_size = 2;
    KmerIndex kmer_index(graph, kmer_size);

    const list<Path> ag_paths = { Path(&graph, 0, { 0 }, 2), Path(&graph, 1, { 2 }, 3) };
    const list<Path> at_paths = { Path(&graph, 0, { 0 }, 2) };

    const list<Path> gg_paths = { Path(&graph, 1, { 0, 1 }, 1), Path(&graph, 0, { 1 }, 2) };
    const list<Path> tg_paths = { Path(&graph, 1, { 0, 1 }, 1) };

    const list<Path> gc_paths = { Path(&graph, 1, { 0, 2 }, 1), Path(&graph, 1, { 1, 2 }, 1) };
    const list<Path> tc_paths = { Path(&graph, 1, { 0, 2 }, 1) };

    const list<Path> ca_paths = { Path(&graph, 0, { 2 }, 2) };

    const StringToPathsMap kmer_to_paths_maps
        = { { "AG", ag_paths }, { "AT", at_paths }, { "GG", gg_paths }, { "TG", tg_paths },
            { "GC", gc_paths }, { "TC", tc_paths }, { "CA", ca_paths } };

    KmerIndex expected_kmer_index(kmer_to_paths_maps);
    ASSERT_EQ(expected_kmer_index, kmer_index);
}

TEST(KmerExtraction, TypicalIndex_KmersExtracted)
{
    Graph graph = makeDeletionGraph("AC", "GG", "CAG");
    const int32_t kmer_size = 2;
    KmerIndex kmer_index(graph, kmer_size);
    const unordered_set<string> expected_kmers = { "AC", "CG", "CC", "GG", "GC", "CA", "AG" };
    ASSERT_EQ(expected_kmers, kmer_index.kmers());
}

TEST(PathExtraction, TypicalIndex_PathsExtracted)
{
    Graph graph = makeDoubleSwapGraph("AAA", "TTT", "CCC", "AAA", "TTT", "AAA", "TTT");
    const int32_t kmer_size = 4;
    KmerIndex kmer_index(graph, kmer_size);
    const list<Path> paths = kmer_index.getPaths("AATT");
    const list<Path> expected_paths
        = { Path(&graph, 1, { 0, 1 }, 2), Path(&graph, 1, { 3, 4 }, 2), Path(&graph, 1, { 5, 6 }, 2) };
    ASSERT_EQ(expected_paths, paths);
}

TEST(CheckingIfKmersArePresent, TypicalKmers_CheckPerformed)
{
    Graph graph = makeDoubleSwapGraph("AAA", "TTT", "CCC", "AAA", "TTT", "AAA", "TTT");
    const int32_t kmer_size = 6;
    KmerIndex kmer_index(graph, kmer_size);
    EXPECT_TRUE(kmer_index.contains("AAATTT"));
    EXPECT_FALSE(kmer_index.contains("AAATTG"));
    EXPECT_FALSE(kmer_index.contains("AAA"));
}

TEST(CountingNumberOfPathsAssociatedWithKmer, TypicalKmers_PathCountObtained)
{
    Graph graph = makeDoubleSwapGraph("AAA", "TTT", "CCC", "AAA", "TTT", "AAA", "TTT");
    {
        const int32_t kmer_size = 6;
        KmerIndex kmer_index(graph, kmer_size);
        EXPECT_EQ(3u, kmer_index.numPaths("AAATTT"));
        EXPECT_EQ(0u, kmer_index.numPaths("AAATTG"));
        EXPECT_EQ(1u, kmer_index.numPaths("TTTTTT"));
    }
    {
        const int32_t kmer_size = 1;
        KmerIndex kmer_index(graph, kmer_size);
        EXPECT_EQ(9u, kmer_index.numPaths("A"));
        EXPECT_EQ(3u, kmer_index.numPaths("C"));
        EXPECT_EQ(9u, kmer_index.numPaths("T"));
        EXPECT_EQ(0u, kmer_index.numPaths("G"));
    }
}

TEST(UniqueKmerCounting, TypicalIndex_UniqueKmersCounted)
{
    Graph graph = makeDeletionGraph("AC", "GG", "ACG");

    const int32_t kmer_size = 3;
    KmerIndex kmer_index(graph, kmer_size);

    EXPECT_EQ(1u, kmer_index.numUniqueKmersOverlappingEdge(0, 1));
    EXPECT_EQ(2u, kmer_index.numUniqueKmersOverlappingEdge(1, 2));

    EXPECT_EQ(3u, kmer_index.numUniqueKmersOverlappingNode(0));
    EXPECT_EQ(4u, kmer_index.numUniqueKmersOverlappingNode(2));
}
