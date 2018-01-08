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

#include "classification/alignment_summary.h"

#include <map>
#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "reads/read.h"

using std::map;
using std::vector;

TEST(SummarizingAlignments, TypicalAlignments_Summarized) {
  reads::ReadPtr flanking_read_1 =
      std::make_shared<reads::Read>("frag1", "AAAACCCCG", "?????????");
  flanking_read_1->SetCanonicalMappingType(MappingType::kFlanksRepeat);
  flanking_read_1->SetNumStrUnitsSpanned(1);

  reads::ReadPtr flanking_read_2 =
      std::make_shared<reads::Read>("frag2", "CCGCCGATT", "?????????");
  flanking_read_2->SetCanonicalMappingType(MappingType::kFlanksRepeat);
  flanking_read_2->SetNumStrUnitsSpanned(2);

  reads::ReadPtr inrepeat_read =
      std::make_shared<reads::Read>("frag3", "CCGCCGCCG", "?????????");
  inrepeat_read->SetCanonicalMappingType(MappingType::kInsideRepeat);
  inrepeat_read->SetNumStrUnitsSpanned(3);

  reads::ReadPtr spanning_read =
      std::make_shared<reads::Read>("frag4", "ACCCCGATT", "?????????");
  spanning_read->SetCanonicalMappingType(MappingType::kSpansRepeat);
  spanning_read->SetNumStrUnitsSpanned(1);

  reads::ReadPtr non_repeat_read =
      std::make_shared<reads::Read>("frag4", "ACTGTGACT", "?????????");
  non_repeat_read->SetCanonicalMappingType(MappingType::kOutsideRepeat);
  non_repeat_read->SetNumStrUnitsSpanned(0);

  vector<reads::ReadPtr> read_ptrs = {flanking_read_1, flanking_read_2,
                                      inrepeat_read, spanning_read,
                                      non_repeat_read};

  map<int32_t, int32_t> flanking_size_counts;
  map<int32_t, int32_t> spanning_size_counts;
  SummarizeAlignments(read_ptrs, flanking_size_counts, spanning_size_counts);

  const map<int32_t, int32_t> expected_flanking_size_counts = {
      {1, 1}, {2, 1}, {3, 1}};
  const map<int32_t, int32_t> expected_spanning_size_counts = {{1, 1}};
  EXPECT_EQ(expected_flanking_size_counts, flanking_size_counts);
  EXPECT_EQ(expected_spanning_size_counts, spanning_size_counts);
}
