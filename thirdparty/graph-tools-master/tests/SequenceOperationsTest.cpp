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
