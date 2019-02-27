//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "common/CountTable.hh"

#include "gtest/gtest.h"

using std::map;
using std::vector;

using namespace ehunter;

TEST(InitializationOfCountTable, TypicalCountTable_Initialized)
{
    const map<int32_t, int32_t> elements_and_counts = { { 1, 2 }, { 3, 5 } };
    CountTable count_table(elements_and_counts);
    EXPECT_EQ(2, count_table.countOf(1));
    EXPECT_EQ(0, count_table.countOf(2));
    EXPECT_EQ(5, count_table.countOf(3));
}

TEST(ManipulatingCountTable, TypicalOperations_TableUpdated)
{
    CountTable count_table;
    count_table.incrementCountOf(4);
    EXPECT_EQ(1, count_table.countOf(4));

    count_table.setCountOf(4, 3);
    EXPECT_EQ(3, count_table.countOf(4));
}

TEST(ObtainingElementsWithNonzeroCounts, TypicalCountTable_ElementsObtained)
{
    const map<int32_t, int32_t> elements_and_counts = { { 1, 2 }, { 3, 5 }, { 7, 15 } };
    CountTable count_table(elements_and_counts);

    count_table.setCountOf(3, 0);

    vector<int32_t> expected_elements = { 1, 7 };
    EXPECT_EQ(expected_elements, count_table.getElementsWithNonzeroCounts());
}
