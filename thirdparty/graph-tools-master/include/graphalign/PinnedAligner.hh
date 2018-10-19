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
