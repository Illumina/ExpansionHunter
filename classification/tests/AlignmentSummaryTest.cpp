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

#include <map>
#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "reads/Read.hh"

using namespace ehunter;

using reads::Read;
using std::map;
using std::vector;

/*
TEST(SummarizingAlignments, TypicalAlignments_Summarized)
{
    Read::SharedPtr flanking_read_1 = std::make_shared<Read>("frag1", "AAAACCCCG");
    flanking_read_1->SetCanonicalAlignmentType(AlignmentType::kFlanksRepeat);
    flanking_read_1->SetNumStrUnitsSpanned(1);

    Read::SharedPtr flanking_read_2 = std::make_shared<Read>("frag2", "CCGCCGATT");
    flanking_read_2->SetCanonicalAlignmentType(AlignmentType::kFlanksRepeat);
    flanking_read_2->SetNumStrUnitsSpanned(2);

    Read::SharedPtr inrepeat_read = std::make_shared<Read>("frag3", "CCGCCGCCG");
    inrepeat_read->SetCanonicalAlignmentType(AlignmentType::kInsideRepeat);
    inrepeat_read->SetNumStrUnitsSpanned(3);

    Read::SharedPtr spanning_read = std::make_shared<Read>("frag4", "ACCCCGATT");
    spanning_read->SetCanonicalAlignmentType(AlignmentType::kSpansRepeat);
    spanning_read->SetNumStrUnitsSpanned(1);

    Read::SharedPtr non_repeat_read = std::make_shared<Read>("frag4", "ACTGTGACT");
    non_repeat_read->SetCanonicalAlignmentType(AlignmentType::kOutsideRepeat);
    non_repeat_read->SetNumStrUnitsSpanned(0);

    vector<Read::SharedPtr> read_ptrs
        = { flanking_read_1, flanking_read_2, inrepeat_read, spanning_read, non_repeat_read };

    CountTable counts_of_flanking_reads;
    CountTable counts_of_spanning_reads;
    SummarizeAlignments(read_ptrs, counts_of_flanking_reads, counts_of_spanning_reads);

    const CountTable expected_counts_of_flanking_reads = { { 1, 1 }, { 2, 1 }, { 3, 1 } };
    const CountTable expected_counts_of_spanning_reads = { { 1, 1 } };
    EXPECT_EQ(expected_counts_of_flanking_reads, counts_of_flanking_reads);
    EXPECT_EQ(expected_counts_of_spanning_reads, counts_of_spanning_reads);
}
*/