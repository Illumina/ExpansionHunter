//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
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

#include "classification/InrepeatReadDetection.hh"

#include <vector>

#include "gmock/gmock.h"

using std::string;
using std::vector;
using testing::DoubleNear;

TEST(CalculatingPeriodicityScore, ValidOffsets_Calculated)
{
    const string sequence = "GGCCCCGGCCCC";
    vector<double> expectedPeriodicityScores = { 0.73, 0.40, 0.33, 0.25, 0.57, 1.00 };

    for (int period = 1; period != 7; ++period)
    {
        const double expectedScore = expectedPeriodicityScores[period - 1];
        EXPECT_THAT(calculatePeriodicityScore(period, sequence), DoubleNear(expectedScore, 0.005));
    }
}

TEST(CalculatingPeriodicityScore, InvalidPeriod_ExceptionThrown)
{
    const string sequence = "GGCCCCGGCCCC";
    ASSERT_ANY_THROW(calculatePeriodicityScore(-1, sequence));
    ASSERT_ANY_THROW(calculatePeriodicityScore(0, sequence));
    ASSERT_ANY_THROW(calculatePeriodicityScore(7, sequence));
}

TEST(DeterminingConsensusRepeatUnit, TypicalSequences_Calculated)
{
    {
        const string sequence = "CGGCGGCGG";
        const int period = 3;
        EXPECT_EQ("CGG", extractConsensusRepeatUnit(period, sequence));
    }
    {
        const string sequence = "CGGATTATTATTCGG";
        const int period = 3;
        EXPECT_EQ("ATT", extractConsensusRepeatUnit(period, sequence));
    }
}

TEST(ComputingMinimalUnitUnderShift, TypicalUnit_Computed)
{
    EXPECT_EQ("CGG", computeSmallestRepeatUnitUnderCircularPermutation("GGC"));
}

TEST(ComputingCanonicalRepeatUnit, TypicalUnit_Computed)
{
    EXPECT_EQ("CCG", computeCanonicalRepeatUnit("CGG"));
    EXPECT_EQ("CCG", computeCanonicalRepeatUnit("GCC"));
}

TEST(ExtractingUnitFromRepetitiveSequences, TypicalRepetitiveSequence_Extracted)
{
    EXPECT_TRUE(checkIfSequenceIsRepetitive("CCG", "CGGCGCCGGCGG"));
    EXPECT_TRUE(checkIfSequenceIsRepetitive("GCG", "CGGCGCCGGCGG"));
    EXPECT_FALSE(checkIfSequenceIsRepetitive("GGG", "CGGCGCCGGCGG"));
    EXPECT_TRUE(checkIfSequenceIsRepetitive("AACCCC", "ACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA"));
    EXPECT_TRUE(checkIfSequenceIsRepetitive("C", "ACCCCACCCCCCCCCCC"));
    EXPECT_FALSE(checkIfSequenceIsRepetitive("A", "ACCCCACCCCCCCCCCC"));
    EXPECT_TRUE(checkIfSequenceIsRepetitive("C", "ACCCCAcccccccccc"));
}