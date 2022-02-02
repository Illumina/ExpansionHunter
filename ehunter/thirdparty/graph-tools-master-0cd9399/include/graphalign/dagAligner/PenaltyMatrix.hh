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

#include <array>
#include <cstdint>

#include <algorithm>
#include <iostream>

#include <boost/assert.hpp>

#include "graphalign/dagAligner/Details.hh"

namespace graphalign
{

namespace dagAligner
{

    class FreePenaltyMatrix
    {
        typedef unsigned char Oligo;

    public:
        typedef Oligo QueryChar;
        typedef Oligo TargetChar;

        // test constructor, since this is a free-form penalty matrix, the one
        // that accepts an actual matrix is expected
        FreePenaltyMatrix(const Score match = 2, const Score mismatch = -2)
            : matrix_{
                // //  a   c   g   t   n
                // {match,    mismatch, mismatch, mismatch, 0}, // a
                // {mismatch, match,    mismatch, mismatch, 0}, // c
                // {mismatch, mismatch, match,    mismatch, 0}, // g
                // {mismatch, mismatch, mismatch, match,    0}, // t
                // {0,        0,        0,        0,        0} // n

                // //  a   c   g   t   n
                // {match,    mismatch, mismatch, mismatch, mismatch}, // a
                // {mismatch, match,    mismatch, mismatch, mismatch}, // c
                // {mismatch, mismatch, match,    mismatch, mismatch}, // g
                // {mismatch, mismatch, mismatch, match,    mismatch}, // t
                // {mismatch, mismatch, mismatch, mismatch, match} // n

                // //  a   c   g   t   n
                // { match, mismatch, mismatch, mismatch, match }, // a
                // { mismatch, match, mismatch, mismatch, match }, // c
                // { mismatch, mismatch, match, mismatch, match }, // g
                // { mismatch, mismatch, mismatch, match, match }, // t
                //

                //  a   c   g   t   n
                { match, mismatch, mismatch, mismatch, match }, // a
                { mismatch, match, mismatch, mismatch, match }, // c
                { mismatch, mismatch, match, mismatch, match }, // g
                { mismatch, mismatch, mismatch, match, match }, // t
                { match, match, match, match, match } // n

            }
        {
        }
        //__attribute__ ((noinline))
        Score operator()(QueryChar q, TargetChar t) const
        {
            if (q >= ROWS_)
            {
                throw std::logic_error("Invalid row request from FreePenaltyMatrix");
            }
            if (t >= COLUMNS_)
            {
                throw std::logic_error("Invalid column request from FreePenaltyMatrix");
            }

            return matrix_[q][t];
        }

        bool isMatch(QueryChar q, TargetChar t) const { return (*this)(q, q) == (*this)(q, t); }

        template <typename InputIterator, typename OutputIterator>
        static void translateTarget(InputIterator targetBegin, InputIterator targetEnd, OutputIterator output)
        {
            std::transform(targetBegin, targetEnd, output, &translateOligo);
        }

        template <typename InputIterator, typename OutputIterator>
        static void translateQuery(InputIterator queryBegin, InputIterator queryEnd, OutputIterator output)
        {
            std::transform(queryBegin, queryEnd, output, &translateOligo);
        }

        friend std::ostream& operator<<(std::ostream& os, const FreePenaltyMatrix& matrix)
        {
            os << "FreePenaltyMatrix(\n";
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

        static const Oligo A = 0;
        static const Oligo C = 1;
        static const Oligo G = 2;
        static const Oligo T = 3;
        static const Oligo N = 4;
        static const Oligo TARGET_CHAR_MAX_ = N;

        static const int ROWS_ = 5;
        static const int COLUMNS_ = 5;

    private:
        const Score matrix_[ROWS_][COLUMNS_];
        static const std::size_t OLIGO_MAX_CHAR = 255;
        static const Oligo TRANSLATION_TABLE_[OLIGO_MAX_CHAR + 1];

        static Oligo translateOligo(unsigned char tc) { return TRANSLATION_TABLE_[tc]; }
    };

    template <bool matchQueryN = true, bool matchTargetN = true> class FixedPenaltyMatrix
    {
        typedef unsigned char Oligo;

    public:
        typedef Oligo QueryChar;
        typedef Oligo TargetChar;

        FixedPenaltyMatrix(const Score match = 2, const Score mismatch = -2)
            : match_(match)
            , mismatch_(mismatch)
        {
        }
        //__attribute__ ((noinline))
        Score operator()(QueryChar q, TargetChar t) const { return isMatch(q, t) ? match_ : mismatch_; }

        bool isMatch(QueryChar q, TargetChar t) const
        {
            return q == t || (matchQueryN && N == q) || (matchTargetN && N == t);
        }

        template <typename InputIterator, typename OutputIterator>
        static void translateTarget(InputIterator targetBegin, InputIterator targetEnd, OutputIterator output)
        {
            std::transform(targetBegin, targetEnd, output, &translateOligo);
        }

        template <typename InputIterator, typename OutputIterator>
        static void translateQuery(InputIterator queryBegin, InputIterator queryEnd, OutputIterator output)
        {
            std::transform(queryBegin, queryEnd, output, &translateOligo);
        }

        friend std::ostream& operator<<(std::ostream& os, const FixedPenaltyMatrix& matrix)
        {
            return os << "FixedPenaltyMatrix(" << matrix.match_ << "," << matrix.mismatch_ << ")";
        }
        static const Oligo A = 0;
        static const Oligo C = 1;
        static const Oligo G = 2;
        static const Oligo T = 3;
        static const Oligo N = 4;
        static const Oligo X = 5;
        static const Oligo TARGET_CHAR_MAX_ = 5;

    private:
        const Score match_;
        const Score mismatch_;
        static const std::size_t OLIGO_MAX_CHAR = 255;
        static const Oligo TRANSLATION_TABLE_[OLIGO_MAX_CHAR + 1];
        static Oligo translateOligo(unsigned char tc) { return TRANSLATION_TABLE_[tc]; }
    };
    // clang-format off
    template <bool matchQueryN, bool matchTargetN>
    const typename FixedPenaltyMatrix<matchQueryN, matchTargetN>::Oligo
        FixedPenaltyMatrix<matchQueryN, matchTargetN>::TRANSLATION_TABLE_[OLIGO_MAX_CHAR + 1]
                                                                           = {
        N, N, N, N, N, N, N, N, N, N,
        N, N, N, N, N, N, N, N, N, N,
        N, N, N, N, N, N, N, N, N, N,
        N, N, N, N, N, N, N, N, N, N,
        N, N, N, N, N, N, N, N, N, N,
        N, N, N, N, N, N, N, N, N, N,
        N, N, N, N, N,
        // capitals
        A, N, C, N, N, N, G, N, N, N,
        N, N, N, N, N, N, N, N, N, T,
        N, N, N, X, N, N,
        // rubbish
        N, N, N, N, N, N,
        // lowercase
        A, N, C, N, N, N, G, N, N, N,
        N, N, N, N, N, N, N, N, N, T,
        N, N, N, X, N, N,
        // padding
        N, N, N, N,
        N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
        N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
        N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
        N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
    };
    // clang-format on

} // namespace dagAligner

} // namespace graphalign
