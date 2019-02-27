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
