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

#include "graphalign/PinnedAligner.hh"

#include <algorithm>

#include "graphalign/TracebackRunner.hh"
#include "graphutils/BaseMatching.hh"

using std::string;

namespace graphtools
{
TracebackMatrix PinnedAligner::populateTracebackMatrix(const string& reference, const string& query)
{
    const size_t num_rows = query.length() + 1;
    const size_t num_cols = reference.length() + 1;

    TracebackMatrix matrix(num_rows, num_cols);

    fillTopLeft(matrix);
    fillTopRow(matrix);
    fillLeftColumn(matrix);
    fillBody(reference, query, matrix);

    return matrix;
}

void PinnedAligner::fillTopLeft(TracebackMatrix& matrix)
{
    matrix.setScore(0, 0, 0);
    matrix.setTracebackStep(0, 0, TracebackStep::kStop);
}

void PinnedAligner::fillTopRow(TracebackMatrix& matrix)
{
    for (size_t col_index = 1; col_index != matrix.numCols(); ++col_index)
    {
        matrix.setScore(0, col_index, col_index * gap_score_);
        matrix.setTracebackStep(0, col_index, TracebackStep::kLeft);
    }
}

void PinnedAligner::fillLeftColumn(TracebackMatrix& matrix)
{
    for (size_t row_index = 1; row_index != matrix.numRows(); ++row_index)
    {
        matrix.setScore(row_index, 0, row_index * gap_score_);
        matrix.setTracebackStep(row_index, 0, TracebackStep::kTop);
    }
}

void PinnedAligner::fillBody(const string& reference, const string& query, TracebackMatrix& matrix)
{
    for (size_t row_index = 1; row_index != matrix.numRows(); ++row_index)
    {
        for (size_t col_index = 1; col_index != matrix.numCols(); ++col_index)
        {
            const bool do_bases_match
                = checkIfReferenceBaseMatchesQueryBase(reference[col_index - 1], query[row_index - 1]);
            fillBodyCell(matrix, row_index, col_index, do_bases_match);
        }
    }
}

void PinnedAligner::fillBodyCell(TracebackMatrix& matrix, size_t row_index, size_t col_index, bool do_bases_match)
{
    const int32_t match_mismatch_score = do_bases_match ? match_score_ : mismatch_score_;
    const TracebackStep traceback_step
        = do_bases_match ? TracebackStep::kDiagonalMatch : TracebackStep::kDiagonalMismatch;
    const int32_t no_gap_score = matrix.score(row_index - 1, col_index - 1) + match_mismatch_score;
    matrix.setScore(row_index, col_index, no_gap_score);
    matrix.setTracebackStep(row_index, col_index, traceback_step);

    const int32_t query_gap_score = matrix.score(row_index, col_index - 1) + gap_score_;
    if (query_gap_score > matrix.score(row_index, col_index))
    {
        matrix.setScore(row_index, col_index, query_gap_score);
        matrix.setTracebackStep(row_index, col_index, TracebackStep::kLeft);
    }

    const int32_t reference_gap_score = matrix.score(row_index - 1, col_index) + gap_score_;
    if (reference_gap_score > matrix.score(row_index, col_index))
    {
        matrix.setScore(row_index, col_index, reference_gap_score);
        matrix.setTracebackStep(row_index, col_index, TracebackStep::kTop);
    }
}

Alignment PinnedAligner::prefixAlign(const string& reference, const string& query)
{
    TracebackMatrix matrix = populateTracebackMatrix(reference, query);

    size_t top_row_index, top_col_index;
    matrix.locateTopScoringCell(top_row_index, top_col_index);

    TracebackRunner traceback_runner(matrix);
    Alignment alignment = traceback_runner.runTraceback(top_row_index, top_col_index);

    return alignment;
}

Alignment PinnedAligner::suffixAlign(string reference, string query)
{
    std::reverse(query.begin(), query.end());
    std::reverse(reference.begin(), reference.end());

    Alignment alignment = prefixAlign(reference, query);
    alignment.reverse(reference.length());

    return alignment;
}
}
