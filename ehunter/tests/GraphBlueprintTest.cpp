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

#include "io/GraphBlueprint.hh"

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
