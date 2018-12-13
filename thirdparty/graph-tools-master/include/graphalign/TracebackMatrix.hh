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

#include <cstdint>
#include <string>
#include <vector>

namespace graphtools
{
enum class TracebackStep
{
    kStop,
    kTop,
    kLeft,
    kDiagonalMatch,
    kDiagonalMismatch
};

std::ostream& operator<<(std::ostream& out, const TracebackStep& trace_dir);

struct TracebackMatrixCell
{
    TracebackMatrixCell()
        : direction(TracebackStep::kStop)
        , score(INT32_MIN)
    {
    }
    TracebackMatrixCell(TracebackStep new_direction, int32_t new_score)
        : direction(new_direction)
        , score(new_score)
    {
    }
    bool operator==(const TracebackMatrixCell& other) const
    {
        return direction == other.direction && score == other.score;
    }
    TracebackStep direction;
    int32_t score;
};

std::ostream& operator<<(std::ostream& out, const TracebackMatrixCell& cell);

/**
 * Implementation of traceback matrix abstraction; each cell of the matrix contains traceback direction and alignment
 * score
 */
class TracebackMatrix
{
public:
    TracebackMatrix(size_t num_rows, size_t num_cols)
        : num_rows_(num_rows)
        , num_cols_(num_cols)
    {
        cells_.resize(num_rows_ * num_cols_);
    };
    explicit TracebackMatrix(const std::string& encoding);

    size_t numRows() const { return num_rows_; }
    size_t numCols() const { return num_cols_; }
    int32_t score(size_t row_index, size_t col_index) const;
    void setScore(size_t row_index, size_t col_index, int32_t score);
    TracebackStep tracebackStep(size_t row_index, size_t col_index) const;
    void setTracebackStep(size_t row_index, size_t col_index, TracebackStep direction);

    void locateTopScoringCell(size_t& top_row_index, size_t& top_col_index) const;

    bool operator==(const TracebackMatrix& other) const;

private:
    inline const TracebackMatrixCell& getCell(size_t row_index, size_t col_index) const
    {
        return cells_[row_index * num_cols_ + col_index];
    }
    inline TracebackMatrixCell& getCell(size_t row_index, size_t col_index)
    {
        return cells_[row_index * num_cols_ + col_index];
    }

    size_t num_rows_;
    size_t num_cols_;
    std::vector<TracebackMatrixCell> cells_;
};

std::ostream& operator<<(std::ostream& out, const TracebackMatrix& matrix);
}
