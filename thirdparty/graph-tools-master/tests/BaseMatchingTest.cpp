//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
