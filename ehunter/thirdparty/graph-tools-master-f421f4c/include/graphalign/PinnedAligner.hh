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

#include "graphalign/LinearAlignment.hh"
#include "graphalign/TracebackMatrix.hh"

namespace graphtools
{

/**
 * Performs local alignment of a pair of sequences that starts at the beginning or the end of both sequences
 */
class PinnedAligner
{
public:
    PinnedAligner(int32_t match_score, int32_t mismatch_score, int32_t gap_score)
        : match_score_(match_score)
        , mismatch_score_(mismatch_score)
        , gap_score_(gap_score)
    {
    }
    TracebackMatrix populateTracebackMatrix(const std::string& reference, const std::string& query);
    // Calculates a top-scoring local alignment of a query to the reference that starts at left-most position of both
    // sequences
    Alignment prefixAlign(const std::string& reference, const std::string& query);
    // Calculates a top-scoring local alignment of a query to the reference that starts at right-most position of both
    // sequences
    Alignment suffixAlign(std::string query, std::string reference);

private:
    void fillTopLeft(TracebackMatrix& matrix);
    void fillTopRow(TracebackMatrix& matrix);
    void fillLeftColumn(TracebackMatrix& matrix);
    void fillBody(const std::string& reference, const std::string& query, TracebackMatrix& matrix);
    void fillBodyCell(TracebackMatrix& matrix, size_t row_index, size_t col_index, bool do_bases_match);

    int32_t match_score_;
    int32_t mismatch_score_;
    int32_t gap_score_;
};
}
