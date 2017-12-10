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

#include "graphs/graph_mapping_operations.h"

#include <string>

#include "gtest/gtest.h"

#include "graphs/graph.h"
#include "graphs/graph_builders.h"

using std::string;

TEST(SplitNodeCigar, ExtractsCigarAndNodeId) {
  const string node_cigar = "1[4M5S]";
  string cigar;
  int32_t node_id;
  splitNodeCigar(node_cigar, cigar, node_id);
  EXPECT_EQ(1, node_id);
  EXPECT_EQ("4M5S", cigar);
}

TEST(DecodeGraphMapping, DecodesTypicalGraphMappings) {
  Graph graph = makeDeletionGraph("AAAA", "TTGG", "TTTT");
  const string read = "AAAATTCCC";
  GraphMapping graph_mapping = decodeFromString(0, "0[4M]1[2M3S]", read, graph);

  GraphMapping expected_graph_mapping(
      {0, 1},
      {Mapping(0, "4M", "AAAA", "AAAA"), Mapping(0, "2M3S", "TTCCC", "TTGG")});
  EXPECT_EQ(expected_graph_mapping, graph_mapping);
}