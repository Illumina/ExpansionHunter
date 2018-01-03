//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "graphs/kmer_index.h"

#include <iostream>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "gtest/gtest.h"

#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/path.h"

using std::cerr;
using std::endl;
using std::list;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using namespace testing;

class TinyDeletionGraph : public ::testing::Test {
 public:
  void SetUp() {
    graph_ptr = MakeDeletionGraph(left_flank, deletion, right_flank);
  }

  const string left_flank = "AC";
  const string deletion = "GG";
  const string right_flank = "CAG";
  GraphSharedPtr graph_ptr;
};

class RepetitiveDoubleSwapGraph : public ::testing::Test {
 protected:
  virtual void SetUp() {
    left_flank = "AAA";
    deletion1 = "TTT";
    insertion1 = "CCC";
    mid = "AAA";
    deletion2 = "TTT";
    insertion2 = "AAA";
    right_flank = "TTT";
    graph_ptr = MakeDoubleSwapGraph(left_flank, deletion1, insertion1, mid,
                                    deletion2, insertion2, right_flank);
  }
  string left_flank;
  string deletion1;
  string insertion1;
  string mid;
  string deletion2;
  string insertion2;
  string right_flank;

  GraphSharedPtr graph_ptr;
};

TEST_F(TinyDeletionGraph, InitializeKmerIndexWith1mers) {
  const int32_t kmer_size = 1;
  KmerIndex kmer_index(graph_ptr, kmer_size);

  const list<GraphPath> a_paths = {GraphPath(graph_ptr, 0, {0}, 0),
                                   GraphPath(graph_ptr, 1, {2}, 1)};
  const list<GraphPath> c_paths = {GraphPath(graph_ptr, 1, {0}, 1),
                                   GraphPath(graph_ptr, 0, {2}, 0)};
  const list<GraphPath> g_paths = {GraphPath(graph_ptr, 0, {1}, 0),
                                   GraphPath(graph_ptr, 1, {1}, 1),
                                   GraphPath(graph_ptr, 2, {2}, 2)};

  const StringToPathsMap kmer_to_paths_maps = {
      {"A", a_paths}, {"C", c_paths}, {"G", g_paths}};

  KmerIndex expected_kmer_index(kmer_to_paths_maps);
  ASSERT_EQ(expected_kmer_index, kmer_index);
}

TEST_F(TinyDeletionGraph, InitializeKmerIndexWith2mers) {
  const int32_t kmer_size = 2;
  KmerIndex kmer_index(graph_ptr, kmer_size);

  const list<GraphPath> ac_paths = {GraphPath(graph_ptr, 0, {0}, 1)};
  const list<GraphPath> cg_paths = {GraphPath(graph_ptr, 1, {0, 1}, 0)};
  const list<GraphPath> cc_paths = {GraphPath(graph_ptr, 1, {0, 2}, 0)};
  const list<GraphPath> gg_paths = {GraphPath(graph_ptr, 0, {1}, 1)};
  const list<GraphPath> gc_paths = {GraphPath(graph_ptr, 1, {1, 2}, 0)};
  const list<GraphPath> ca_paths = {GraphPath(graph_ptr, 0, {2}, 1)};
  const list<GraphPath> ag_paths = {GraphPath(graph_ptr, 1, {2}, 2)};

  const StringToPathsMap kmer_to_paths_maps = {
      {"AC", ac_paths}, {"CG", cg_paths}, {"CC", cc_paths}, {"GG", gg_paths},
      {"GC", gc_paths}, {"CA", ca_paths}, {"AG", ag_paths}};

  KmerIndex expected_kmer_index(kmer_to_paths_maps);
  ASSERT_EQ(expected_kmer_index, kmer_index);
}

TEST_F(TinyDeletionGraph, KmerIndexReportsKmersWithNonzeroCount) {
  const int32_t kmer_size = 2;
  KmerIndex kmer_index(graph_ptr, kmer_size);
  const unordered_set<string> expected_kmers = {"AC", "CG", "CC", "GG",
                                                "GC", "CA", "AG"};
  ASSERT_EQ(expected_kmers, kmer_index.GetKmersWithNonzeroCount());
}

TEST_F(RepetitiveDoubleSwapGraph, ExtractPathsContainingKmer) {
  const int32_t kmer_size = 4;
  KmerIndex kmer_index(graph_ptr, kmer_size);
  const list<GraphPath> paths = kmer_index.GetPaths("AATT");
  const list<GraphPath> expected_paths = {GraphPath(graph_ptr, 1, {0, 1}, 1),
                                          GraphPath(graph_ptr, 1, {3, 4}, 1),
                                          GraphPath(graph_ptr, 1, {5, 6}, 1)};
  ASSERT_EQ(expected_paths, paths);
}

TEST_F(RepetitiveDoubleSwapGraph, CheckIfIndexContainsKmer) {
  const int32_t kmer_size = 6;
  KmerIndex kmer_index(graph_ptr, kmer_size);
  EXPECT_TRUE(kmer_index.Contains("AAATTT"));
  EXPECT_FALSE(kmer_index.Contains("AAATTG"));
  EXPECT_FALSE(kmer_index.Contains("AAA"));
}

TEST_F(RepetitiveDoubleSwapGraph, GetNumberOfPathsContainingKmer) {
  {
    const int32_t kmer_size = 6;
    KmerIndex kmer_index(graph_ptr, kmer_size);
    EXPECT_EQ(3u, kmer_index.NumPaths("AAATTT"));
    EXPECT_EQ(0u, kmer_index.NumPaths("AAATTG"));
    EXPECT_EQ(1u, kmer_index.NumPaths("TTTTTT"));
  }
  {
    const int32_t kmer_size = 1;
    KmerIndex kmer_index(graph_ptr, kmer_size);
    EXPECT_EQ(9u, kmer_index.NumPaths("A"));
    EXPECT_EQ(3u, kmer_index.NumPaths("C"));
    EXPECT_EQ(9u, kmer_index.NumPaths("T"));
    EXPECT_EQ(0u, kmer_index.NumPaths("G"));
  }
}