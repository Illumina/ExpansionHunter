//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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

#include <iostream>

#include "graphalign/dagAligner/Details.hh"
#include "graphutils/BaseMatching.hh"

namespace graphalign
{

namespace dagAligner
{

    class BaseMatchingPenaltyMatrix
    {
        typedef graphtools::codes::BaseCode Oligo;
        typedef graphalign::dagAligner::Score Score;

    public:
        typedef Oligo QueryChar;
        typedef Oligo TargetChar;

        // test constructor, since this is a free-form penalty matrix, the one
        // that accepts an actual matrix is expected
        BaseMatchingPenaltyMatrix(const Score match = 2, const Score mismatch = -2)
        {
            int row = 0;
            for (const auto& r : graphtools::codes::kReferenceQueryCodeMatchLookupTable)
            {
                int col = 0;
                for (const auto& c : r)
                {
                    matrix_[row][col] = c ? match : mismatch;
                    ++col;
                }
                ++row;
            }
        }

        Score operator()(QueryChar q, TargetChar t) const
        {
            if (t >= ROWS_)
            {
                throw std::logic_error("Invalid row request from BaseMatchingPenaltyMatrix");
            }
            if (q >= COLUMNS_)
            {
                throw std::logic_error("Invalid column request from BaseMatchingPenaltyMatrix");
            }

            return matrix_[t][q];
        }

        bool isMatch(QueryChar q, TargetChar t) const
        {
            return graphtools::codes::kReferenceQueryCodeMatchLookupTable[t][q];
        }

        template <typename InputIterator, typename OutputIterator>
        static void translateTarget(InputIterator targetBegin, InputIterator targetEnd, OutputIterator output)
        {
            std::transform(targetBegin, targetEnd, output, [](unsigned char tc) {
                return graphtools::codes::kReferenceBaseEncodingTable[tc];
            });
        }

        template <typename InputIterator, typename OutputIterator>
        static void translateQuery(InputIterator queryBegin, InputIterator queryEnd, OutputIterator output)
        {
            std::transform(queryBegin, queryEnd, output, [](unsigned char tc) {
                return graphtools::codes::kQueryBaseEncodingTable[tc];
            });
        }

        friend std::ostream& operator<<(std::ostream& os, const BaseMatchingPenaltyMatrix& matrix)
        {
            os << "BaseMatchingPenaltyMatrix(\n";
            for (int r = 0; r < ROWS_; ++r)
            {
                for (int c = 0; c < COLUMNS_; ++c)
                {
                    os << (c ? "\t" : "[") << int(matrix.matrix_[r][c]);
                }
                os << "]\n";
            }
            return os << ")";
        }

        static const int ROWS_ = graphtools::codes::kMaxReferenceBaseCode + 1;
        static const int COLUMNS_ = graphtools::codes::kMaxQueryBaseCode + 1;
        static const Oligo TARGET_CHAR_MAX_ = graphtools::codes::kMaxReferenceBaseCode;

    private:
        graphalign::dagAligner::Score matrix_[ROWS_][COLUMNS_];
    };

} // namespace dagAligner

} // namespace graphalign
