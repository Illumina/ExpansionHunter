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

#include "classification/overlap_quantification.h"

#include <string>

#include "gtest/gtest.h"

#include "graphs/graph.h"
#include "graphs/graph_builders.h"
#include "graphs/graph_mapping.h"
#include "graphs/graph_mapping_operations.h"

using std::string;

TEST(StrOverlapQuantification, TypicalReads_StrOverlapComputed) {
  GraphSharedPtr graph_ptr = MakeStrGraph("ATAT", "CCG", "ATTT");
  int32_t str_unit_len = 3;
  StrOverlapQuantifier str_overlap_quantifier(0, 1, 2, str_unit_len);

  {
    const string non_repeat_read = "ATAT";
    GraphMapping mapping =
        DecodeFromString(0, "0[4M]", non_repeat_read, graph_ptr);
    const int32_t num_units =
        str_overlap_quantifier.NumUnitsOverlapped(mapping);
    ASSERT_EQ(0, num_units);
  }

  {
    const string spanning_read = "ATCCGCCGAT";
    GraphMapping mapping =
        DecodeFromString(2, "0[2M]1[3M]1[3M]2[2M]", spanning_read, graph_ptr);
    const int32_t num_units =
        str_overlap_quantifier.NumUnitsOverlapped(mapping);
    ASSERT_EQ(2, num_units);
  }

  {
    const string flanking_read = "ATCCGCCGCC";
    GraphMapping mapping =
        DecodeFromString(2, "0[2M]1[3M]1[3M]1[2M]", flanking_read, graph_ptr);
    const int32_t num_units =
        str_overlap_quantifier.NumUnitsOverlapped(mapping);
    ASSERT_EQ(2, num_units);
  }

  {
    const string inrepeat_read = "CCGCCGCCGCC";
    GraphMapping mapping =
        DecodeFromString(0, "1[3M]1[3M]1[3M]1[2M]", inrepeat_read, graph_ptr);
    const int32_t num_units =
        str_overlap_quantifier.NumUnitsOverlapped(mapping);
    ASSERT_EQ(3, num_units);
  }
}
