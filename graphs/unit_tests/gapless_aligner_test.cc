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

#include "graphs/gapless_aligner.h"

#include <memory>
#include <string>

#include "gtest/gtest.h"

#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/graph_mapping.h"
#include "graphs/graph_mapping_operations.h"
#include "graphs/path.h"

using std::string;

class DeletionGraphForAlignment : public ::testing::Test {
 protected:
  virtual void SetUp() {
    left_flank = "AAAACC";
    deletion = "TTTGG";
    right_flank = "ATTT";
    graph = makeDeletionGraph(left_flank, deletion, right_flank);
    graph_ptr = std::make_shared<Graph>(graph);
  }
  Graph graph;
  string left_flank;
  string deletion;
  string right_flank;
  std::shared_ptr<Graph> graph_ptr;
};

TEST_F(DeletionGraphForAlignment, AligningReadsToShortReferenceCausesError) {
  EXPECT_ANY_THROW(alignWithoutGaps("AAAA", 0, "AAA"));
}

TEST_F(DeletionGraphForAlignment, AligningEmptySequencesCausesError) {
  EXPECT_ANY_THROW(alignWithoutGaps("", 0, ""));
}

TEST_F(DeletionGraphForAlignment, AlignSequences) {
  const string query = "AGGTTTTG";
  const string reference = "NNNNATCGTTTG";
  const Mapping expected_mapping(4, "1M3X4M", query, reference);
  ASSERT_EQ(expected_mapping, alignWithoutGaps(query, 4, reference));
}

TEST_F(DeletionGraphForAlignment, AlignReadToSingleNodePath) {
  GraphPath path(graph_ptr, 1, {1}, 4);
  const string read = "ATGC";

  GraphMapping expected_graph_mapping =
      decodeFromString(1, "1[1X2M1X]", read, graph);
  GraphMapping graph_mapping = alignWithoutGaps(path, read);
  EXPECT_EQ(expected_graph_mapping, graph_mapping);
}

TEST_F(DeletionGraphForAlignment, AlignReadToMultiNodePath) {
  GraphPath path(graph_ptr, 2, {0, 1, 2}, 1);
  const string read = "TTCCTTAGGAT";

  GraphMapping expected_graph_mapping =
      decodeFromString(2, "0[2X2M]1[2M1X2M]2[2M]", read, graph);
  GraphMapping graph_mapping = alignWithoutGaps(path, read);
  EXPECT_EQ(expected_graph_mapping, graph_mapping);
}