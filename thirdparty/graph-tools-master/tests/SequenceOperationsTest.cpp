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

#include "graphutils/SequenceOperations.hh"

#include <map>
#include <string>
#include <vector>

#include "gtest/gtest.h"

using std::map;
using std::string;
using std::vector;

using namespace graphtools;

TEST(CheckingSequenceComposition, TypicalSequences_CompositionDetermined)
{
    const string reference_nucleotide_sequence = "ACTG";
    const string reference_sequence = "ACWG";
    const string nonreference_sequence = "ZZZZ";

    EXPECT_TRUE(checkIfNucleotideReferenceSequence(reference_nucleotide_sequence));
    EXPECT_FALSE(checkIfNucleotideReferenceSequence(reference_sequence));

    EXPECT_TRUE(checkIfReferenceSequence(reference_nucleotide_sequence));
    EXPECT_TRUE(checkIfReferenceSequence(reference_sequence));
    EXPECT_FALSE(checkIfReferenceSequence(nonreference_sequence));
}

TEST(ExpandingDegenerateSymbols, TypicalSymbol_SymbolExpanded)
{
    map<char, string> kSymbolExpansion
        = { { 'A', "A" },   { 'C', "C" },   { 'T', "T" },    { 'G', "G" },  { 'R', "AG" },  { 'Y', "CT" },
            { 'K', "GT" },  { 'M', "AC" },  { 'S', "CG" },   { 'W', "AT" }, { 'B', "CGT" }, { 'D', "AGT" },
            { 'H', "ACT" }, { 'V', "ACG" }, { 'N', "ACGT" }, { 'X', "X" } };

    for (const auto& symbol_expansion : kSymbolExpansion)
    {
        EXPECT_EQ(symbol_expansion.second, expandReferenceSymbol(symbol_expansion.first));
    }
}

TEST(ExpandingDegenerateSymbols, NonReferenceSymbol_ExceptionThrown) { ASSERT_ANY_THROW(expandReferenceSymbol('a')); }

TEST(ExpandingDegenerateSequences, SequenceWithDegenerateBases_SequenceExpanded)
{
    string sequence = "RAK";
    const vector<string> expected_expansion = { "AAG", "GAG", "AAT", "GAT" };
    vector<string> observed_expansion = { "AAG", "GAG", "AAT", "GAT" };
    expandReferenceSequence(sequence, observed_expansion);
    ASSERT_EQ(expected_expansion, observed_expansion);
}

TEST(SplittingStrings, WordsDelimitedBySpaces_StringVector)
{
    const string composite_string = "abc /+=  ##";
    const vector<string> expected_words = { "abc", "/+=", "##" };
    ASSERT_EQ(expected_words, splitStringByWhitespace(composite_string));
}

TEST(SplittingStrings, WordsDelimitedBySlashes_StringVector)
{
    const string string_with_words = "a/b/cd";
    const vector<string> expected_words = { "a", "b", "cd" };
    ASSERT_EQ(expected_words, splitStringByDelimiter(string_with_words, '/'));
}

TEST(ReverseComplementingSequences, TypicalQueryAndReferenceSequences_ReverseComplemented)
{
    EXPECT_EQ("AAGGCGAT", reverseComplement("ATCGCCTT"));
    EXPECT_EQ("aaggcgat", reverseComplement("atcgcctt"));
    EXPECT_EQ("RYKMSWBDHVN", reverseComplement("NBDHVWSKMRY"));
}
