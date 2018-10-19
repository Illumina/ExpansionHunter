// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

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