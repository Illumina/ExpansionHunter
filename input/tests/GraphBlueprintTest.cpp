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

#include "input/GraphBlueprint.hh"

#include "gtest/gtest.h"

using std::string;
using std::vector;

using namespace ehunter;

TEST(SplittingStringsIntoTokens, ValidStrings_Split)
{
    const string regex = "ATGC(CAG)+GTCG(AAA|TTT)(AGTC)?(CAG)*";
    auto tokens = tokenizeRegex(regex);

    vector<string> expectedTokens = { "ATGC", "(CAG)+", "GTCG", "(AAA|TTT)", "(AGTC)?", "(CAG)*" };
    ASSERT_EQ(expectedTokens, tokens);
}

TEST(ParsingTokens, TypicalTokens_Parsed)
{
    TokenParser parser;
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kInsertionOrDeletion, { "AGTC" } };
        EXPECT_EQ(expectedResult, parser.parse("(AGTC)?"));
    }
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kSkippableRepeat, { "CAG" } };
        EXPECT_EQ(expectedResult, parser.parse("(CAG)*"));
    }
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kUnskippableRepeat, { "CAG" } };
        EXPECT_EQ(expectedResult, parser.parse("(CAG)+"));
    }
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kInterruption, { "GTCG" } };
        EXPECT_EQ(expectedResult, parser.parse("GTCG"));
    }
    {
        FeatureTypeAndSequences expectedResult = { GraphBlueprintFeatureType::kSwap, { "AAA", "TTT" } };
        EXPECT_EQ(expectedResult, parser.parse("(AAA|TTT)"));
    }
}