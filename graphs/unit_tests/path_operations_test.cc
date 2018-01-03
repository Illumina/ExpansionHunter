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

#include "graphs/path_operations.h"

#include "gtest/gtest.h"

#include <list>
#include <string>
#include <vector>

#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/path.h"

using std::list;
using std::string;
using std::vector;

TEST(SplittingSequenceByPath, SequenceOfDifferentLength_ExceptionRaised) {
  GraphSharedPtr graph_ptr = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");
  GraphPath path(graph_ptr, 3, {0, 1}, 2);
  const string sequence = "AA";
  EXPECT_ANY_THROW(SplitByPath(path, sequence));
}

TEST(SplittingSequenceByPath, SingleNodePath_SequenceSplit) {
  GraphSharedPtr graph_ptr = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");
  GraphPath path(graph_ptr, 1, {1}, 3);
  const string sequence = "AAT";
  const vector<string> expected_pieces = {sequence};
  EXPECT_EQ(expected_pieces, SplitByPath(path, sequence));
}

TEST(SplittingSequenceByPath, MultiNodePath_SequenceSplit) {
  GraphSharedPtr graph_ptr = MakeDeletionGraph("AAAACC", "TTTGG", "ATTT");
  {
    GraphPath path(graph_ptr, 1, {0, 1}, 3);
    const string sequence = "AAAAAGGGG";
    const vector<string> expected_pieces = {"AAAAA", "GGGG"};
    EXPECT_EQ(expected_pieces, SplitByPath(path, sequence));
  }
  {
    GraphPath path(graph_ptr, 3, {0, 2}, 1);
    const string sequence = "AAACC";
    const vector<string> expected_pieces = {"AAA", "CC"};
    EXPECT_EQ(expected_pieces, SplitByPath(path, sequence));
  }
  {
    GraphPath path(graph_ptr, 3, {0, 1, 2}, 1);
    const string sequence = "AAAGGGGGCC";
    const vector<string> expected_pieces = {"AAA", "GGGGG", "CC"};
    EXPECT_EQ(expected_pieces, SplitByPath(path, sequence));
  }
}