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

#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace graphtools
{

using StringPair = std::pair<std::string, std::string>;

/**
 * Splits a string by the specified delimiter
 */
std::vector<std::string> splitStringByDelimiter(const std::string& str, char sep = ' ');

/**
 * Splits a string by whitespace
 */
std::vector<std::string> splitStringByWhitespace(const std::string& str);

/**
 * Return reversed sequence
 * @param seq sequence to reverse
 * @return reversed sequence
 */
static inline std::string reverseString(std::string seq)
{
    std::reverse(seq.begin(), seq.end());
    return seq;
}

/**
 * Return reverse complemented sequence
 * @param seq sequence to reverse-complement
 * @return reverse complemented sequence
 */
std::string reverseComplement(std::string seq);

/**
 * Returns true if sequence consists of uppercase symbols over extended nucleotide alphabet
 */
bool checkIfReferenceSequence(const std::string& sequence);

/**
 * Checks if sequence consists of uppercase As, Ts, Cs, and Gs
 */
bool checkIfNucleotideReferenceSequence(const std::string& sequence);

/**
 * Expands reference symbol into strings made up of matching nucleotides
 */
std::string const& expandReferenceSymbol(char reference_symbol);

/**
 * Expands reference sequence by expanding each degenerate symbol
 *
 * @param sequence Reference sequence to expand
 * @return All nucleotide sequences matching the input sequence
 */
void expandReferenceSequence(const std::string& sequence, std::vector<std::string>& target);
}
