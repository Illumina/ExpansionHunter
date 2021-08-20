//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Roman Petrovski <RPetrovski@illumina.com>
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
