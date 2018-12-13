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

using namespace ehunter;

TEST(SearchingForHighQualityBaseRuns, AllBasesHighQuality_FullRangeReturned)
{
    string sequence = "ATCGATCG";
    const auto goodBasesRange = findHighQualityBaseRun(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin(), sequence.end());
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, SequenceEndingInLowQualityBases_CorrectRangeReturned)
{
    string sequence = "ATCGATCgaTcg";
    const auto goodBasesRange = findHighQualityBaseRun(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin(), sequence.end() - 5);
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, SequenceStartingWithLowQualityBases_CorrectRangeReturned)
{
    string sequence = "gaTcgaTCGATC";
    const auto goodBasesRange = findHighQualityBaseRun(sequence, 0.1, 0.8);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin() + 5, sequence.end());
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, SequenceFlankedByLowQualityBasesOnBothSides_CorrectRangeReturned)
{
    string sequence = "gaTcgaTCGATCgaTcg";
    const auto goodBasesRange = findHighQualityBaseRun(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin() + 6, sequence.end() - 5);
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, SequenceComprisedOfLowQualityBases_EmptyRangeReturned)
{
    string sequence = "gaTcgaatgtTCatg";
    const auto goodBasesRange = findHighQualityBaseRun(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin() + 6, sequence.end() - 3);
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, RealReadEndingInManyLowQualityBases_CorrectRangeReturned)
{
    const string sequence = "CCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACaGCCGCCACCGCCGCCGCCGCCGCC"
                            "GCCGCCtCCgCAGCCtCCtCaGCCGCCGCCGCCgcCgCaGCCGCcGCcgCCgCcgcCgcc";

    const auto goodBasesRange = findHighQualityBaseRun(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin(), sequence.end() - 7);
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, QueryWithLowQualityStart_CorrectRangeReturned)
{
    const string sequence = "GcgggggGcGgcggcggcGggggcgcgggggccgGggggcGtGCGGcgggggggcGGcGGcGGCGGggGCGGcGGcGGcGGCGGcGgCGG"
                            "CGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG";

    const auto goodBasesRange = findHighQualityBaseRun(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin() + 55, sequence.end());
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, QueryWithLowQualityEnd_CorrectRangeReturned)
{
    const string sequence = "GGcGGcGGCGGggGCGGcGGcGGcGGCGGcGgCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGC"
                            "GGCGGGcgggggGcGgcggcggcGggggcgcgggggccgGggggcGtGCGGcgggggggc";

    const auto goodBasesRange = findHighQualityBaseRun(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin(), sequence.end() - 54);
    ASSERT_EQ(expectedBases, goodBases);
}

TEST(SearchingForHighQualityBaseRuns, QueryWithLowQualityEnds_CorrectRangeReturned)
{
    const string sequence
        = "cgggggccgGggggcGtGCGGcgggggGGcGGcGGCGGggGCGGcGGcGGcGGCGGcGgCGGCGGCGGCGGCGGCGGCGGCGGGGGCGGGA"
          "cgggggGcGgcggcggcGggggcgcgggggccgGggggcGtGCGGcgggggggc";

    const auto goodBasesRange = findHighQualityBaseRun(sequence);
    string goodBases = string(goodBasesRange.first, goodBasesRange.second);

    string expectedBases = string(sequence.begin() + 27, sequence.end() - 54);
    ASSERT_EQ(expectedBases, goodBases);
}
