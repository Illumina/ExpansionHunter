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

TEST(AligningSequences, SequencesWithUnequalLength_ExceptionThrown) {
  EXPECT_ANY_THROW(AlignWithoutGaps("AAAA", 0, "AAA"));
}

TEST(AligningSequences, EmptySequences_ExceptionThrown) {
  EXPECT_ANY_THROW(AlignWithoutGaps("", 0, ""));
}

TEST(AligningSequences, TypicalSequences_Aligned) {
  const string query = "AGGTTTTG";
  const string reference = "NNNNATCGTTTG";
  const Mapping expected_mapping(4, "1M3X4M", query, reference);
  ASSERT_EQ(expected_mapping, AlignWithoutGaps(query, 4, reference));
}

TEST(AligningSequenceToPath, SingleNodePath_Aligned) {
  Graph graph = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  GraphPath path(graph_ptr, 1, {1}, 4);
  const string read = "ATGC";

  GraphMapping expected_graph_mapping =
      DecodeFromString(1, "1[1X2M1X]", read, graph);
  GraphMapping graph_mapping = AlignWithoutGaps(path, read);
  EXPECT_EQ(expected_graph_mapping, graph_mapping);
}

TEST(AligningSequenceToPath, MultiNodePath_Aligned) {
  Graph graph = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");
  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
  GraphPath path(graph_ptr, 2, {0, 1, 2}, 1);
  const string read = "TTCCTTAGGAT";

  GraphMapping expected_graph_mapping =
      DecodeFromString(2, "0[2X2M]1[2M1X2M]2[2M]", read, graph);
  GraphMapping graph_mapping = AlignWithoutGaps(path, read);
  EXPECT_EQ(expected_graph_mapping, graph_mapping);
}

// TEST(AligningSequenceToPath, TypicalStrPath_Aligned) {
//  Graph graph = MakeStrGraph("AAAACC", "CCG", "ATTT");
//  std::shared_ptr<Graph> graph_ptr = std::make_shared<Graph>(graph);
//  GraphPath path(graph_ptr, 2, {0, 1, 2}, 1);
//  //                   FFFFRRRRRRRRRFFFF
//  const string read = "AACCCCGCCGCCGATTT";
//
//  GraphMapping expected_graph_mapping =
//      DecodeFromString(2, "0[4M]1[3M]1[3M]1[3M]2[4M]", read, graph);
//  GraphMapping graph_mapping = AlignWithoutGaps(path, read);
//  EXPECT_EQ(expected_graph_mapping, graph_mapping);
//}