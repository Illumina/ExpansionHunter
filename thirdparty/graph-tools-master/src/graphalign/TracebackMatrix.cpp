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

#include "graphalign/TracebackMatrix.hh"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "graphutils/SequenceOperations.hh"

using std::string;
using std::vector;

namespace graphtools
{

TracebackStep decodeTracebackDirection(std::string encoding)
{
    if (encoding == "S")
    {
        return TracebackStep::kStop;
    }
    else if (encoding == "I")
    {
        return TracebackStep::kTop;
    }
    else if (encoding == "D")
    {
        return TracebackStep::kLeft;
    }
    else if (encoding == "M")
    {
        return TracebackStep::kDiagonalMatch;
    }
    else if (encoding == "X")
    {
        return TracebackStep::kDiagonalMismatch;
    }
    throw std::logic_error(encoding + " is an unknown traceback");
}

TracebackMatrix::TracebackMatrix(const string& encoding)
{
    vector<string> lines = splitStringByDelimiter(encoding, '\n');
    num_rows_ = lines.size();
    num_cols_ = std::count(lines.front().begin(), lines.front().end(), '/');
    cells_.resize(num_rows_ * num_cols_);

    for (size_t row_index = 0; row_index != num_rows_; ++row_index)
    {
        const std::string& line = lines[row_index];
        std::vector<std::string> words = splitStringByWhitespace(line);
        for (size_t col_index = 0; col_index != num_cols_; ++col_index)
        {
            const std::string& word = words[col_index];
            std::vector<std::string> tracedir_score_encodings = splitStringByDelimiter(word, '/');
            assert(tracedir_score_encodings.size() == 2);
            const TracebackStep trace_dir = decodeTracebackDirection(tracedir_score_encodings.front());
            const int32_t score = std::stoi(tracedir_score_encodings.back());
            getCell(row_index, col_index) = TracebackMatrixCell(trace_dir, score);
        }
    }
}

int32_t TracebackMatrix::score(size_t row_index, size_t col_index) const { return getCell(row_index, col_index).score; }

void TracebackMatrix::setScore(size_t row_index, size_t col_index, int32_t score)
{
    getCell(row_index, col_index).score = score;
}

TracebackStep TracebackMatrix::tracebackStep(size_t row_index, size_t col_index) const
{
    return getCell(row_index, col_index).direction;
}

void TracebackMatrix::setTracebackStep(size_t row_index, size_t col_index, TracebackStep direction)
{
    getCell(row_index, col_index).direction = direction;
}

void TracebackMatrix::locateTopScoringCell(size_t& top_row_index, size_t& top_col_index) const
{
    int32_t top_score = INT32_MIN;

    for (size_t row_index = 0; row_index != num_rows_; ++row_index)
    {
        for (size_t col_index = 0; col_index != num_cols_; ++col_index)
        {
            const int32_t cell_score = score(row_index, col_index);
            if (top_score <= cell_score)
            {
                top_score = cell_score;
                top_row_index = row_index;
                top_col_index = col_index;
            }
        }
    }
}

bool TracebackMatrix::operator==(const TracebackMatrix& other) const
{
    return (num_rows_ == other.num_rows_ && num_cols_ == other.num_cols_ && cells_ == other.cells_);
}

std::ostream& operator<<(std::ostream& out, const TracebackStep& trace_dir)
{
    switch (trace_dir)
    {
    case TracebackStep::kStop:
        out << "S";
        break;
    case TracebackStep::kTop:
        out << "I";
        break;
    case TracebackStep::kLeft:
        out << "D";
        break;
    case TracebackStep::kDiagonalMatch:
        out << "M";
        break;
    case TracebackStep::kDiagonalMismatch:
        out << "X";
        break;
    default:
        out << "?";
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const TracebackMatrixCell& cell)
{
    out << cell.direction << "/" << cell.score;
    return out;
}

std::ostream& operator<<(std::ostream& out, const TracebackMatrix& matrix)
{
    out << "{";
    for (size_t row_index = 0; row_index != matrix.numRows(); ++row_index)
    {
        out << "{";
        for (size_t col_index = 0; col_index != matrix.numCols(); ++col_index)
        {
            out << matrix.tracebackStep(row_index, col_index) << "/" << matrix.score(row_index, col_index);
            if (col_index != matrix.numCols() - 1)
            {
                out << ", ";
            }
        }
        if (row_index != matrix.numRows() - 1)
        {
            out << "}, ";
        }
    }
    out << "}}";
    return out;
}
}
