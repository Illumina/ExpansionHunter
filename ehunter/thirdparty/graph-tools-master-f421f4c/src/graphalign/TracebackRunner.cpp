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
