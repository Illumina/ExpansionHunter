//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "common/count_table.h"

#include "gtest/gtest.h"

using std::map;
using std::vector;

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
