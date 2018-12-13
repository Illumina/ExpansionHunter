//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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

#include "graphalign/dagAligner/PenaltyMatrix.hh"

namespace graphalign
{

namespace dagAligner
{

    // clang-format off
// Notice: this one does not support X at the moment
const FreePenaltyMatrix::Oligo FreePenaltyMatrix::TRANSLATION_TABLE_[OLIGO_MAX_CHAR + 1] = {
    N, N, N, N, N, N, N, N, N, N,
    N, N, N, N, N, N, N, N, N, N,
    N, N, N, N, N, N, N, N, N, N,
    N, N, N, N, N, N, N, N, N, N,
    N, N, N, N, N, N, N, N, N, N,
    N, N, N, N, N, N, N, N, N, N,
    N, N, N, N, N,
    // capitals
    A, N, C, N, N, N, G, N, N, N,
    N, N, N, N, N, N, N, N, N, T,
    N, N, N, N, N, N,
    // rubbish
    N, N, N, N, N, N,
    // lowercase
    A, N, C, N, N, N, G, N, N, N,
    N, N, N, N, N, N, N, N, N, T,
    N, N, N, N, N, N,
    // padding
    N, N, N, N,
    N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
    N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
    N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
    N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N, N,
};
    // clang-format on

} // namespace dagAligner

} // namespace graphalign
