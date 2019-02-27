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
