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

#include "graphalign/TracebackMatrix.hh"

#include <string>

#include "gtest/gtest.h"

using std::string;

using namespace graphtools;

TEST(TracebackMatrixInitialization, DefaultInitialization_UnirectedMatrixOfZeros)
{
    TracebackMatrix traceback_matrix(2, 3);

    for (size_t row_index = 0; row_index != 2; ++row_index)
    {
        for (size_t col_index = 0; col_index != 3; ++col_index)
        {
            EXPECT_EQ(INT32_MIN, traceback_matrix.score(row_index, col_index));
            EXPECT_EQ(TracebackStep::kStop, traceback_matrix.tracebackStep(row_index, col_index));
        }
    }
}

TEST(TracebackMatrixInitialization, TypicalEncoding_MatrixInitialized)
{
    const std::string encoding = "S/0  D/-2  D/-4\n"
                                 "I/0  M/-1  I/-4";

    TracebackMatrix traceback_matrix(encoding);

    TracebackMatrix expected_matrix(2, 3);
    expected_matrix.setScore(0, 0, 0);
    expected_matrix.setTracebackStep(0, 0, TracebackStep::kStop);

    expected_matrix.setScore(0, 1, -2);
    expected_matrix.setTracebackStep(0, 1, TracebackStep::kLeft);

    expected_matrix.setScore(0, 2, -4);
    expected_matrix.setTracebackStep(0, 2, TracebackStep::kLeft);

    expected_matrix.setScore(1, 0, 0);
    expected_matrix.setTracebackStep(1, 0, TracebackStep::kTop);

    expected_matrix.setScore(1, 1, -1);
    expected_matrix.setTracebackStep(1, 1, TracebackStep::kDiagonalMatch);

    expected_matrix.setScore(1, 2, -4);
    expected_matrix.setTracebackStep(1, 2, TracebackStep::kTop);

    ASSERT_EQ(expected_matrix, traceback_matrix);
}

TEST(LocatingTopScoringCell, TypicalMatrix_CellLocated)
{
    TracebackMatrix matrix("S/0   S/2  S/10\n"
                           "S/-1  S/3 S/-1");

    size_t row_index;
    size_t col_index;
    matrix.locateTopScoringCell(row_index, col_index);
    EXPECT_EQ(0u, row_index);
    EXPECT_EQ(2u, col_index);
}
