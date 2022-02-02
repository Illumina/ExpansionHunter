//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>
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

    const std::vector<Path> expected_a_paths = { Path(&graph, 0, { 0 }, 1), Path(&graph, 1, { 2 }, 2) };
    const std::vector<Path> expected_c_paths = { Path(&graph, 1, { 0 }, 2), Path(&graph, 0, { 2 }, 1) };
    const std::vector<Path> expected_g_paths
        = { Path(&graph, 0, { 1 }, 1), Path(&graph, 1, { 1 }, 2), Path(&graph, 2, { 2 }, 3) };

    ASSERT_EQ(kmer_index.getPaths("A"), expected_a_paths);
    ASSERT_EQ(kmer_index.getPaths("C"), expected_c_paths);
    ASSERT_EQ(kmer_index.getPaths("G"), expected_g_paths);
}

TEST(KmerIndexInitialization, 2mers_IndexCreated)
{
    Graph graph = makeDeletionGraph("AK", "GG", "CAG");
    const int32_t kmer_size = 2;
    KmerIndex kmer_index(graph, kmer_size);

    const std::vector<Path> expected_ag_paths = { Path(&graph, 0, { 0 }, 2), Path(&graph, 1, { 2 }, 3) };
    const std::vector<Path> expected_at_paths = { Path(&graph, 0, { 0 }, 2) };

    const std::vector<Path> expected_gg_paths = { Path(&graph, 1, { 0, 1 }, 1), Path(&graph, 0, { 1 }, 2) };
    const std::vector<Path> expected_tg_paths = { Path(&graph, 1, { 0, 1 }, 1) };

    const std::vector<Path> expected_gc_paths = { Path(&graph, 1, { 0, 2 }, 1), Path(&graph, 1, { 1, 2 }, 1) };
    const std::vector<Path> expected_tc_paths = { Path(&graph, 1, { 0, 2 }, 1) };

    const std::vector<Path> expected_ca_paths = { Path(&graph, 0, { 2 }, 2) };

    ASSERT_EQ(kmer_index.getPaths("AG"), expected_ag_paths);
    ASSERT_EQ(kmer_index.getPaths("AT"), expected_at_paths);
    ASSERT_EQ(kmer_index.getPaths("GG"), expected_gg_paths);
    ASSERT_EQ(kmer_index.getPaths("TG"), expected_tg_paths);
    ASSERT_EQ(kmer_index.getPaths("GC"), expected_gc_paths);
    ASSERT_EQ(kmer_index.getPaths("TC"), expected_tc_paths);
    ASSERT_EQ(kmer_index.getPaths("CA"), expected_ca_paths);
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
    const std::vector<Path> paths = kmer_index.getPaths("AATT");
    const std::vector<Path> expected_paths
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
