//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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
