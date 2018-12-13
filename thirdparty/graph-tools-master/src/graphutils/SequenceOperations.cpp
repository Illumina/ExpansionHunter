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

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

using std::string;
using std::unordered_map;
using std::vector;

namespace graphtools
{

vector<string> splitStringByDelimiter(const std::string& str, char sep)
{
    vector<string> words;
    std::stringstream sstream(str);
    string word;
    while (std::getline(sstream, word, sep))
    {
        words.push_back(word);
    }

    return words;
}

vector<string> splitStringByWhitespace(const std::string& str)
{
    vector<string> words;
    std::stringstream sstream(str);
    string word;
    while (sstream >> word)
    {
        words.push_back(word);
    }
    return words;
}

static char complementBase(char base)
{
    switch (base)
    {
    case 'A':
        return 'T';
    case 'a':
        return 't';
    case 'C':
        return 'G';
    case 'c':
        return 'g';
    case 'G':
        return 'C';
    case 'g':
        return 'c';
    case 'T':
        return 'A';
    case 't':
        return 'a';
    case 'R':
        return 'Y';
    case 'Y':
        return 'R';
    case 'K':
        return 'M';
    case 'M':
        return 'K';
    case 'S':
        return 'S';
    case 'W':
        return 'W';
    case 'B':
        return 'V';
    case 'D':
        return 'H';
    case 'H':
        return 'D';
    case 'V':
        return 'B';
    default:
        return 'N';
    }
}

string reverseComplement(string seq)
{
    std::transform(seq.begin(), seq.end(), seq.begin(), complementBase);
    std::reverse(seq.begin(), seq.end());
    return seq;
}

const unordered_map<char, string> kSymbolExpansion
    = { { 'A', "A" },   { 'C', "C" },   { 'T', "T" },    { 'G', "G" },  { 'R', "AG" },  { 'Y', "CT" },
        { 'K', "GT" },  { 'M', "AC" },  { 'S', "CG" },   { 'W', "AT" }, { 'B', "CGT" }, { 'D', "AGT" },
        { 'H', "ACT" }, { 'V', "ACG" }, { 'N', "ACGT" }, { 'X', "X" } };

static bool checkIfNucleotideReferenceSymbol(char symbol)
{
    return (symbol == 'A') || (symbol == 'C') || (symbol == 'T') || (symbol == 'G');
}

static bool hasExpandableSymbols(const std::string& s)
{
    static const struct ExpandableSyms
    {
        ExpandableSyms()
        {
            for (const auto& sym : kSymbolExpansion)
            {
                if (sym.second.size() > 1)
                {
                    value.push_back(sym.first);
                }
            }
        }
        string value;
    } expandableSyms;
    return s.find_first_of(expandableSyms.value) != string::npos;
}

bool checkIfNucleotideReferenceSequence(const std::string& sequence)
{
    for (char symbol : sequence)
    {
        if (!checkIfNucleotideReferenceSymbol(symbol))
        {
            return false;
        }
    }
    return true;
}

static bool checkIfReferenceSymbol(char symbol) { return kSymbolExpansion.find(symbol) != kSymbolExpansion.end(); }

bool checkIfReferenceSequence(const std::string& sequence)
{
    for (char symbol : sequence)
    {
        if (!checkIfReferenceSymbol(symbol))
        {
            return false;
        }
    }
    return true;
}

string const& expandReferenceSymbol(char symbol)
{
    if (!checkIfReferenceSymbol(symbol))
    {
        const string symbol_encoding(1, symbol);
        throw std::logic_error("Symbol " + symbol_encoding + " is not a valid reference symbol");
    }
    return kSymbolExpansion.at(symbol);
}

void expandReferenceSequence(const string& sequence, vector<string>& expanded_sequences)
{
    if (!hasExpandableSymbols(sequence))
    {
        expanded_sequences = { sequence };
        return;
    }
    expanded_sequences.resize(1);
    expanded_sequences.front().clear();
    expanded_sequences.front().reserve(2 * sequence.size());

    for (char symbol : sequence)
    {
        const auto& expansions = expandReferenceSymbol(symbol);

        // first expansion: append to all
        for (auto& expanded_sequence : expanded_sequences)
        {
            expanded_sequence.push_back(expansions.front());
        }

        size_t sequences_to_expand = expanded_sequences.size();
        size_t exp_pos = 1;
        while (exp_pos < expansions.size())
        {
            // more than one expansion: create expanded copy for each of them
            for (size_t to_copy_pos = 0; to_copy_pos < sequences_to_expand; ++to_copy_pos)
            {
                string to_copy = expanded_sequences[to_copy_pos];
                to_copy.back() = expansions[exp_pos];
                expanded_sequences.emplace_back(std::move(to_copy));
            }
            exp_pos++;
        }
    }
}
}
