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
