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

#include "graphalign/TracebackRunner.hh"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>

using std::string;
using std::vector;

namespace graphtools
{

Alignment TracebackRunner::runTraceback(size_t row_index, size_t col_index)
{
    operations_.clear();

    if (row_index != matrix_.numRows() - 1)
    {
        softclipQuerySuffix(row_index);
    }

    while (matrix_.tracebackStep(row_index, col_index) != TracebackStep::kStop)
    {
        computeTracebackRunForAlignmentOperation(row_index, col_index);
        convertCurrentRunToAlignmentOperation();
        row_index = run_last_row_index;
        col_index = run_last_col_index;
        tracebackPosition(row_index, col_index);
    }

    if (row_index != 0)
    {
        softclipQueryPrefix(row_index);
    }

    std::reverse(operations_.begin(), operations_.end());
    return Alignment(col_index, operations_);
}

void TracebackRunner::computeTracebackRunForAlignmentOperation(size_t row_index, size_t col_index)
{
    TracebackStep cur_traceback_step = matrix_.tracebackStep(row_index, col_index);
    run_traceback_step = cur_traceback_step;
    run_length = 0;

    while (run_traceback_step == cur_traceback_step)
    {
        run_last_col_index = col_index;
        run_last_row_index = row_index;
        ++run_length;
        tracebackPosition(row_index, col_index);
        cur_traceback_step = matrix_.tracebackStep(row_index, col_index);
    }
}

void TracebackRunner::tracebackPosition(size_t& row_index, size_t& col_index)
{
    TracebackStep traceback_step = matrix_.tracebackStep(row_index, col_index);
    switch (traceback_step)
    {
    case TracebackStep::kDiagonalMatch:
    case TracebackStep::kDiagonalMismatch:
        --row_index;
        --col_index;
        break;
    case TracebackStep::kLeft:
        --col_index;
        break;
    case TracebackStep::kTop:
        --row_index;
        break;
    case TracebackStep::kStop:
        break;
    }
}

void TracebackRunner::convertCurrentRunToAlignmentOperation()
{
    OperationType operation_type;

    switch (run_traceback_step)
    {
    case TracebackStep::kDiagonalMatch:
        operation_type = OperationType::kMatch;
        break;
    case TracebackStep::kDiagonalMismatch:
        operation_type = OperationType::kMismatch;
        break;
    case TracebackStep::kLeft:
        operation_type = OperationType::kDeletionFromRef;
        break;
    case TracebackStep::kTop:
        operation_type = OperationType::kInsertionToRef;
        break;
    case TracebackStep::kStop:
    default:
        throw std::logic_error("Attempted invalid traceback run decoding");
    }

    operations_.emplace_back(operation_type, run_length);
}

void TracebackRunner::softclipQuerySuffix(size_t& row_index)
{

    const uint32_t softclip_len = matrix_.numRows() - row_index - 1;

    operations_.emplace_back(OperationType::kSoftclip, softclip_len);
}

void TracebackRunner::softclipQueryPrefix(size_t& row_index)
{
    operations_.emplace_back(OperationType::kSoftclip, row_index);
}
}
