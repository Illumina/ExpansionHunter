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

#include "alignment/HighQualityBaseRunFinder.hh"

#include <string>

#include "gtest/gtest.h"

using std::make_pair;
using std::string;

TEST(SearchingForHighQualityBaseRuns, AllBasesHighQuality_FullRangeReturned)
{
    const int windowSize = 6;
    const int minHighQualityBasesInGoodWindow = 3;
    const int minLengthOfFullSizeRun = 5;
    HighQualityBaseRunFinder goodBaseFinder(windowSize, minHighQualityBasesInGoodWindow, minLengthOfFullSizeRun);

    string sequence = "ATCGATCG";
    const auto goodBasesRange = goodBaseFinder.find(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin(), sequence.end());
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, SequenceEndingInLowQualityBases_CorrectRangeReturned)
{
    HighQualityBaseRunFinder goodBaseFinder(6, 3, 5);

    string sequence = "ATCGATCgaTcg";
    const auto goodBasesRange = goodBaseFinder.find(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin(), sequence.end() - 6);
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, SequenceStartingWithLowQualityBases_CorrectRangeReturned)
{
    HighQualityBaseRunFinder goodBaseFinder(6, 3, 5);

    string sequence = "gaTcgaTCGATC";
    const auto goodBasesRange = goodBaseFinder.find(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin() + 2, sequence.end());
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, SequenceFlankedByLowQualityBasesOnBothSides_CorrectRangeReturned)
{
    HighQualityBaseRunFinder goodBaseFinder(6, 3, 5);

    string sequence = "gaTcgaTCGATCgaTcg";
    const auto goodBasesRange = goodBaseFinder.find(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin() + 2, sequence.end() - 6);
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, SequenceComprisedOfLowQualityBases_EmptyRangeReturned)
{
    HighQualityBaseRunFinder goodBaseFinder(6, 3, 5);

    string sequence = "gaTcgaatgtTCatg";
    const auto goodBasesRange = goodBaseFinder.find(sequence);

    HighQualityBaseRunFinder::StringIterPair expectedRange = make_pair(sequence.end(), sequence.end());
    ASSERT_EQ(expectedRange, goodBasesRange);
}

TEST(SearchingForHighQualityBaseRuns, RealReadEndingInManyLowQualityBases_CorectRangeReturned)
{
    HighQualityBaseRunFinder goodBaseFinder(6, 3, 75);

    const string sequence = "CCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACaGCCGCCACCGCCGCCGCCGCCGCC"
                            "GCCGCCtCCgCAGCCtCCtCaGCCGCCGCCGCCgcCgCaGCCGCcGCcgCCgCcgcCgcc";
    const auto goodBasesRange = goodBaseFinder.find(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin(), sequence.end() - 27);
    ASSERT_EQ(expectedBases, goodBases);
}
