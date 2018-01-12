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

#include "common/seq_operations.h"
#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/graph_mapping.h"
#include "graphs/graph_mapping_operations.h"
#include "graphs/path.h"

using std::list;
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
  GraphSharedPtr graph_ptr = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");
  GraphPath path(graph_ptr, 1, {1}, 4);
  const string read = "ATGC";

  GraphMapping expected_graph_mapping =
      DecodeFromString(1, "1[1X2M1X]", read, graph_ptr);
  GraphMapping graph_mapping = AlignWithoutGaps(path, read);
  EXPECT_EQ(expected_graph_mapping, graph_mapping);
}

TEST(AligningSequenceToPath, MultiNodePath_Aligned) {
  GraphSharedPtr graph_ptr = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");
  GraphPath path(graph_ptr, 2, {0, 1, 2}, 1);
  const string read = "TTCCTTAGGAT";

  GraphMapping expected_graph_mapping =
      DecodeFromString(2, "0[2X2M]1[2M1X2M]2[2M]", read, graph_ptr);
  GraphMapping graph_mapping = AlignWithoutGaps(path, read);
  EXPECT_EQ(expected_graph_mapping, graph_mapping);
}

TEST(AligningSequenceToPath, TypicalStrPath_Aligned) {
  GraphSharedPtr graph_ptr = MakeStrGraph("AAAACC", "CCG", "ATTT");
  GraphPath path(graph_ptr, 2, {0, 1, 1, 1, 2}, 3);
  //                   FFFFRRRRRRRRRFFFF
  const string read = "AACCCCGCCGCCGATTT";

  GraphMapping expected_graph_mapping =
      DecodeFromString(2, "0[4M]1[3M]1[3M]1[3M]2[4M]", read, graph_ptr);
  GraphMapping graph_mapping = AlignWithoutGaps(path, read);
  EXPECT_EQ(expected_graph_mapping, graph_mapping);
}

TEST(KmerExtraction, TypicalSequence_KmersExtracted) {
  const string sequence = "AAATTT";
  const list<string> expected_4mers = {"AAAT", "AATT", "ATTT"};
  ASSERT_EQ(expected_4mers, ExtractKmersFromAllPositions(sequence, 4));

  const list<string> expected_7mers = {};
  ASSERT_EQ(expected_7mers, ExtractKmersFromAllPositions(sequence, 7));
}

TEST(AlignmentOfSequenceToShortPath, TypicalSequence_BestAlignmentObtained) {
  GraphSharedPtr graph_ptr = MakeDeletionGraph("AAACC", "TTGGG", "TTAAA");
  const GraphPath path(graph_ptr, 4, {0}, 4);
  const string sequence = "CCTTA";

  list<GraphMapping> mappings = GetBestAlignmentToShortPath(path, 1, sequence);

  list<GraphMapping> expected_mappings = {
      DecodeFromString(3, "0[2M]2[3M]", sequence, graph_ptr)};
  ASSERT_EQ(expected_mappings, mappings);
}

TEST(AlignmentOfSequenceToGraph, TypicalSequence_BestAlignmentObtained) {
  GraphSharedPtr graph_ptr = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");

  const int32_t kmer_len = 3;
  GaplessAligner aligner(graph_ptr, kmer_len);
  const string sequence = "TTCCTTAGGAT";

  list<GraphMapping> mappings = aligner.GetBestAlignment(sequence);

  list<GraphMapping> expected_mappings = {
      DecodeFromString(2, "0[2X2M]1[2M1X2M]2[2M]", sequence, graph_ptr)};
  ASSERT_EQ(expected_mappings, mappings);
}

TEST(GraphAlignment, TypicalStrGraph_BestAlignmentObtained) {
  GraphSharedPtr graph_ptr = MakeStrGraph("AAAACG", "CCG", "ATTT");
  const int32_t kmer_len = 3;
  GaplessAligner aligner(graph_ptr, kmer_len);

  {
    //                            FFFFRRRRRRRRRFFFF
    const string spanning_read = "AACGCCGCCGCCGATTT";
    list<GraphMapping> mappings = aligner.GetBestAlignment(spanning_read);

    list<GraphMapping> expected_mappings = {DecodeFromString(
        2, "0[4M]1[3M]1[3M]1[3M]2[4M]", spanning_read, graph_ptr)};
    EXPECT_EQ(expected_mappings, mappings);
  }

  {
    //                          RRRRRRRRRRR
    const string repeat_read = "CGCCGCCGCCG";
    list<GraphMapping> mappings = aligner.GetBestAlignment(repeat_read);
    list<GraphMapping> expected_mappings = {
        DecodeFromString(4, "0[2M]1[3M]1[3M]1[3M]", repeat_read, graph_ptr),
        DecodeFromString(1, "1[2M]1[3M]1[3M]1[3M]", repeat_read, graph_ptr)};
    EXPECT_EQ(expected_mappings, mappings);
  }

  {
    //                          RRRXRRRRXRRR
    const string repeat_read = "CCGACGCCTCCG";
    list<GraphMapping> mappings = aligner.GetBestAlignment(repeat_read);
    list<GraphMapping> expected_mappings = {DecodeFromString(
        0, "1[3M]1[1X2M]1[2M1X]1[3M]", repeat_read, graph_ptr)};
    EXPECT_EQ(expected_mappings, mappings);
  }
}

TEST(StrandClassification, TypicalRead_StrandDetermined) {
  GraphSharedPtr graph_ptr = MakeStrGraph("AAAACC", "CCG", "ATTT");
  const int32_t kmer_len = 3;
  StrandClassifier classifier(graph_ptr, kmer_len);
  EXPECT_TRUE(classifier.IsForwardOriented("CCGCCGCCGCCG"));
  EXPECT_FALSE(classifier.IsForwardOriented(ReverseComplement("CCGCCGCCGCCG")));
  EXPECT_TRUE(classifier.IsForwardOriented("CCGACGCCTCCG"));
  EXPECT_FALSE(classifier.IsForwardOriented(ReverseComplement("CCGACGCCTCCG")));
}
