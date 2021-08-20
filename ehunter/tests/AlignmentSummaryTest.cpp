//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

#include <map>
#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "core/Read.hh"

using namespace ehunter;

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
