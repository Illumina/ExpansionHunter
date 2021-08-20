//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
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

#include "graphutils/BaseMatching.hh"

#include "gtest/gtest.h"

using namespace graphtools;

TEST(MatchingBases, CoreBases_Matched)
{
    EXPECT_TRUE(checkIfReferenceBaseMatchesQueryBase('A', 'A'));
    EXPECT_TRUE(checkIfReferenceBaseMatchesQueryBase('T', 't'));
    EXPECT_TRUE(checkIfReferenceBaseMatchesQueryBase('c', 'C'));
    EXPECT_TRUE(checkIfReferenceBaseMatchesQueryBase('G', 'g'));

    EXPECT_FALSE(checkIfReferenceBaseMatchesQueryBase('T', 'c'));
}

TEST(MatchingBases, DegenerateBases_Matched)
{
    EXPECT_TRUE(checkIfReferenceBaseMatchesQueryBase('Y', 'c'));
    EXPECT_TRUE(checkIfReferenceBaseMatchesQueryBase('Y', 'T'));
    EXPECT_TRUE(checkIfReferenceBaseMatchesQueryBase('W', 'a'));
    EXPECT_TRUE(checkIfReferenceBaseMatchesQueryBase('N', 'A'));

    EXPECT_FALSE(checkIfReferenceBaseMatchesQueryBase('Y', 'a'));
}

TEST(MatchingStrings, TypicalStrings_Matched)
{
    EXPECT_TRUE(checkIfReferenceAndQuerySequencesMatch("RRTCS", "AaTCG"));
    EXPECT_FALSE(checkIfReferenceAndQuerySequencesMatch("WG", "CG"));
    EXPECT_FALSE(checkIfReferenceAndQuerySequencesMatch("CC", "CCC"));
}
