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

#include "reads/read_operations.h"

#include "gtest/gtest.h"

#include "graphs/gapless_aligner.h"
#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "reads/read.h"

using namespace reads;

TEST(ReorientingReads, CorrectlyOrientedRead_ReadIsUnchnaged) {
  Read read("frag1", "AAATCT", "(((???");

  Read expected_read = read;

  Graph graph = MakeStrGraph("AAAA", "CG", "TCTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  const int32_t kmer_len = 3;
  StrandClassifier classifier(graph_ptr, kmer_len);

  ReorientRead(classifier, read);

  ASSERT_EQ(expected_read, read);
}

TEST(ReorientingReads, IncorrectlyOrientedRead_BasesAndQualsReoriented) {
  Read read("frag1", "GACGTT", "?#?(((");

  Graph graph = MakeStrGraph("AAAA", "CG", "TCTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  const int32_t kmer_len = 3;
  StrandClassifier classifier(graph_ptr, kmer_len);

  ReorientRead(classifier, read);

  Read expected_read("frag1", "AACGTC", "(((?#?");

  ASSERT_EQ(expected_read, read);
}
