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

#pragma once

#include <array>
#include <cassert>
#include <string>

namespace graphtools
{

namespace codes
{
    using BaseCode = unsigned char;
    const int kMaxBaseCode = 15;

    // Core base codes
    const BaseCode A = 0;
    const BaseCode C = 1;
    const BaseCode G = 2;
    const BaseCode T = 3;
    const BaseCode X = 4;

    // Degenerate base codes
    const BaseCode B = 5;
    const BaseCode D = 6;
    const BaseCode H = 7;
    const BaseCode K = 8;
    const BaseCode M = 9;
    const BaseCode N = 10;
    const BaseCode R = 11;
    const BaseCode S = 12;
    const BaseCode V = 13;
    const BaseCode W = 14;
    const BaseCode Y = 15;

    const int kMaxQueryBaseCode = 4;
    const int kMaxReferenceBaseCode = 15;

    const int maxBaseAscii = 255;

    // Currently, low-quality (lower case) bases get the same encoding as their high-quality counterparts. We should
    // extend the coding scheme when we are ready to deal with base-quality in the alignment.

    // Core bases A, C, G, T and degenerate bases B, D, H, K, M, N, S, R, V, W, Y all receive distinct codes. All other
    // base symbols are coded as X, which is the code intended to mismatch everything.
    const std::array<BaseCode, maxBaseAscii + 1> kReferenceBaseEncodingTable
        = { X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, A, B, C, D, X, X, G, H, X, X, K, X, M, N, X, X, X, R, S, T, X, V, W, X, Y, X, X, X, X, X, X,
            X, A, X, C, X, X, X, G, X, X, X, X, X, X, X, X, X, X, X, X, T, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X };

    // Core bases A, C, G, T all recieve distinct codes. All other base symbols are coded as X.
    const std::array<BaseCode, maxBaseAscii + 1> kQueryBaseEncodingTable
        = { X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, A, X, C, X, X, X, G, X, X, X, X, X, X, X, X, X, X, X, X, T, X, X, X, X, X, X, X, X, X, X, X,
            X, A, X, C, X, X, X, G, X, X, X, X, X, X, X, X, X, X, X, X, T, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,
            X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X };

    // We use the standard matching rules for degenerate bases. The X symbol corresponds to a mismatch
    const std::array<std::array<bool, codes::kMaxQueryBaseCode + 1>, codes::kMaxReferenceBaseCode + 1>
        kReferenceQueryCodeMatchLookupTable = {
            // clang-format off
          //   A      C      G      T      X
          { { true,  false, false, false, false }, // A
            { false, true,  false, false, false }, // C
            { false, false, true,  false, false }, // G
            { false, false, false, true,  false }, // T
            { false, false, false, false, false }, // X
            { false, true,  true,  true,  false }, // B
            { true,  false, true,  true,  false }, // D
            { true,  true,  false, true,  false }, // H
            { false, false, true,  true,  false }, // K
            { true,  true,  false, false, false }, // M
            { true,  true,  true,  true,  false }, // N
            { true,  false, true,  false, false }, // R
            { false, true,  true,  false, false }, // S
            { true,  true,  true,  false, false }, // V
            { true,  false, false, true,  false }, // W
            { false, true,  false, true,  false }  // Y
            }
            // clang-format on
        };
}

inline codes::BaseCode encodeReferenceBase(char base) { return codes::kReferenceBaseEncodingTable[base]; }

inline codes::BaseCode encodeQueryBase(char base) { return codes::kQueryBaseEncodingTable[base]; }

/**
 * Checks if a pair of reference and query base codes corresponds to matching bases
 *
 * Examples:
 *   checkIfReferenceBaseCodeMatchesQueryBaseCode(encodeReferenceBase('C'), encodeQueryBase('c')) == true
 *   checkIfReferenceBaseCodeMatchesQueryBaseCode(encodeReferenceBase('T'), encodeQueryBase('Y')) == true
 *   checkIfReferenceBaseCodeMatchesQueryBaseCode(encodeReferenceBase('a'), encodeQueryBase('W')) == true
 *   checkIfReferenceBaseCodeMatchesQueryBaseCode(encodeReferenceBase('C'), encodeQueryBase('G')) == false
 *
 */
inline bool checkIfReferenceBaseCodeMatchesQueryBaseCode(codes::BaseCode referenceCode, codes::BaseCode queryCode)
{
    // Temporary asserts until we add validation code for reference/query strings.
    assert(referenceCode <= codes::kMaxReferenceBaseCode);
    assert(queryCode <= codes::kMaxQueryBaseCode);

    return codes::kReferenceQueryCodeMatchLookupTable[referenceCode][queryCode];
}

inline bool checkIfReferenceBaseMatchesQueryBase(char referenceBase, char queryBase)
{
    return checkIfReferenceBaseCodeMatchesQueryBaseCode(encodeReferenceBase(referenceBase), encodeQueryBase(queryBase));
}

/**
 * Checks if a reference sequence matches a query sequence
 *
 * @param reference And reference sequence
 * @param query Any query sequence
 * @return True if the sequences match
 */
inline bool checkIfReferenceAndQuerySequencesMatch(const std::string& reference, const std::string& query)
{
    if (reference.length() != query.length())
    {
        return false;
    }

    for (size_t index = 0; index != reference.length(); ++index)
    {
        if (!checkIfReferenceBaseMatchesQueryBase(reference[index], query[index]))
        {
            return false;
        }
    }

    return true;
}
}
