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

#include "graphs/graph_mapping.h"

#include "gtest/gtest.h"

#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/graph_mapping_operations.h"

using std::list;
using std::string;
using std::vector;

TEST(Operation, InitializesFromString) {
  Operation operation("3M", "ATC", "ATC");
  Operation expected_operation('M', 3, "ATC", "ATC");
  ASSERT_EQ(expected_operation, operation);
}

TEST(Operation, OutputsQueryAndReferenceSpans) {
  {
    Operation operation("3M", "AAA", "AAA");
    EXPECT_EQ((int32_t)3, operation.QuerySpan());
    EXPECT_EQ((int32_t)3, operation.ReferenceSpan());
  }
  {
    Operation operation("4X", "AAAA", "TTTT");
    EXPECT_EQ((int32_t)4, operation.QuerySpan());
    EXPECT_EQ((int32_t)4, operation.ReferenceSpan());
  }
  {
    Operation operation("5D", "", "AAAAA");
    EXPECT_EQ((int32_t)0, operation.QuerySpan());
    EXPECT_EQ((int32_t)5, operation.ReferenceSpan());
  }
  {
    Operation operation("7I", "AAAAAAA", "");
    EXPECT_EQ((int32_t)7, operation.QuerySpan());
    EXPECT_EQ((int32_t)0, operation.ReferenceSpan());
  }
  {
    Operation operation("10S", "AAAAAAAAAA", "");
    EXPECT_EQ((int32_t)10, operation.QuerySpan());
    EXPECT_EQ((int32_t)0, operation.ReferenceSpan());
  }
  {
    Operation operation("7N", "NNNNNNN", "NNNNNNN");
    EXPECT_EQ((int32_t)7, operation.QuerySpan());
    EXPECT_EQ((int32_t)7, operation.ReferenceSpan());
  }
  {
    Operation operation("3N", "NCN", "CNN");
    EXPECT_EQ((int32_t)3, operation.QuerySpan());
    EXPECT_EQ((int32_t)3, operation.ReferenceSpan());
  }
}

TEST(Operation, ErrorsOutOnUnexpectedSequences) {
  EXPECT_ANY_THROW(Operation("4M", "AAAA", "ATCG"));
  EXPECT_ANY_THROW(Operation("4M", "AAAA", "ATC"));
  EXPECT_ANY_THROW(Operation("4M", "AAA", "AAA"));

  EXPECT_ANY_THROW(Operation("4N", "NNN", "NNN"));
  EXPECT_ANY_THROW(Operation("3N", "NN", "NNN"));
  EXPECT_ANY_THROW(Operation("2N", "NT", "NT"));

  EXPECT_ANY_THROW(Operation("2X", "AT", "TT"));
  EXPECT_ANY_THROW(Operation("2X", "AT", "A"));

  EXPECT_ANY_THROW(Operation("4D", "AAA", ""));
  EXPECT_ANY_THROW(Operation("4D", "", ""));

  EXPECT_ANY_THROW(Operation("2I", "AA", "T"));

  EXPECT_ANY_THROW(Operation("2S", "TTT", ""));
  EXPECT_ANY_THROW(Operation("2S", "TT", "T"));
}

TEST(Mapping, InitializesFromCigar) {
  // query: ---TTCGTT--TTGGGTCCCCCCCCCC
  //           ||| ||  ||   |
  //   ref: CCCTTCCNNAATT---T----------
  string query = "TTCGTTTTGGGTCCCCCCCCCC";
  string reference = "CCCTTCCNNAATTT";

  Mapping mapping(3, "3M1X2N2D2M3I1M10S", query, reference);
  const vector<Operation> operations = {
      Operation('M', 3, "TTC", "TTC"), Operation('X', 1, "G", "C"),
      Operation('N', 2, "TT", "NN"),   Operation('D', 2, "", "AA"),
      Operation('M', 2, "TT", "TT"),   Operation('I', 3, "GGG", ""),
      Operation('M', 1, "T", "T"),     Operation('S', 10, "CCCCCCCCCC", "")};

  Mapping expected_mapping(3, operations);
  ASSERT_EQ(expected_mapping, mapping);
}

TEST(Mapping, CalculatesQueryAndReferenceSpans) {
  Mapping mapping(3, "3M1X2M2D2M3I1M10S", "TTCGTTTTGGGTCCCCCCCCCC",
                  "CCCTTCCTTAATTT");
  EXPECT_EQ((int32_t)22, mapping.QuerySpan());
  EXPECT_EQ((int32_t)11, mapping.ReferenceSpan());
}

TEST(GraphMapping, CalculatesNumberOfMatches) {
  Graph graph = MakeDeletionGraph("AAAA", "TTGG", "TTTT");
  const string query = "AAAATTCCC";
  GraphMapping graph_mapping =
      DecodeFromString(0, "0[4M]1[2M3S]", query, graph);
  EXPECT_EQ((int32_t)6, graph_mapping.NumMatches());
}

TEST(Mapping, OutputsQueryAndReferenceSequences) {
  Mapping mapping(3, "3M1X2M2D2M3I1M10S", "TTCGTTTTGGGTCCCCCCCCCC",
                  "CCCTTCCTTAATTT");
  EXPECT_EQ("TTCGTTTTGGGT", mapping.Query());
  EXPECT_EQ("TTCCTTAATTT", mapping.Reference());
}

TEST(GraphMapping, StitchesQueryAndReferenceSequences) {
  Graph graph = MakeDeletionGraph("AAAA", "TTGG", "TTTT");
  const string query = "AAAATTCCC";
  GraphMapping graph_mapping =
      DecodeFromString(0, "0[4M]1[2M3S]", query, graph);
  EXPECT_EQ("AAAATT", graph_mapping.Query());
  EXPECT_EQ("AAAATT", graph_mapping.Reference());
}

TEST(GraphMapping, CalculatesQueryAndReferenceSpans) {
  Graph graph = MakeDeletionGraph("AAAA", "TTGG", "TTTT");
  const string query = "AAAATTCCC";
  GraphMapping graph_mapping =
      DecodeFromString(0, "0[4M]1[2M3S]", query, graph);
  EXPECT_EQ((int32_t)9, graph_mapping.QuerySpan());
  EXPECT_EQ((int32_t)6, graph_mapping.ReferenceSpan());
}

TEST(GraphMapping, AllowsAccessingNodeMappingsByIndex) {
  Graph graph = MakeDeletionGraph("AAAA", "TTGC", "TTTT");
  const string query = "AAAATTCCC";
  GraphMapping graph_mapping =
      DecodeFromString(0, "0[4M]1[2M3S]", query, graph);
  EXPECT_EQ(Mapping(0, "4M", "AAAA", "AAAA"), graph_mapping[0].mapping);
  EXPECT_EQ(Mapping(0, "2M3S", "TTCCC", "TTGG"), graph_mapping[1].mapping);
}

TEST(GettingIndexesOfNode, TypicalMapping_IndexesObtained) {
  Graph graph = MakeStrGraph("AAAACC", "CCG", "ATTT");
  const string read = "CCCCGCCGAT";
  GraphMapping mapping =
      DecodeFromString(4, "0[2M]1[3M]1[3M]2[2M]", read, graph);
  const list<int32_t> left_flank_indexes = {0};
  const list<int32_t> repeat_unit_indexes = {1, 2};
  const list<int32_t> right_flank_indexes = {3};
  EXPECT_EQ(left_flank_indexes, mapping.GetIndexesOfNode(0));
  EXPECT_EQ(repeat_unit_indexes, mapping.GetIndexesOfNode(1));
  EXPECT_EQ(right_flank_indexes, mapping.GetIndexesOfNode(2));
}

TEST(GettingIndexesOfNode, NodeNotInMapping_EmptyListReturned) {
  Graph graph = MakeStrGraph("AAAACC", "CCG", "ATTT");
  const string read = "ACCCCG";
  GraphMapping mapping = DecodeFromString(3, "0[3M]1[3M]", read, graph);
  const list<int32_t> empty_list;
  EXPECT_EQ(empty_list, mapping.GetIndexesOfNode(2));
  EXPECT_EQ(empty_list, mapping.GetIndexesOfNode(4));
}

TEST(CheckingIfMappingSpansNode, TypicalMapping_ChecksPerformed) {
  Graph graph = MakeStrGraph("AAAACC", "CCG", "ATTT");
  const string read = "ACCCCG";
  GraphMapping mapping = DecodeFromString(3, "0[3M]1[3M]", read, graph);
  EXPECT_TRUE(mapping.SpansNode(0));
  EXPECT_TRUE(mapping.SpansNode(1));
  EXPECT_FALSE(mapping.SpansNode(2));
  EXPECT_FALSE(mapping.SpansNode(3));
}