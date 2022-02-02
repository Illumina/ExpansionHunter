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

#include <cassert>
#include <cstdint>
/**
 * Holds scores for linear sequence alignment algorithms
 */
struct LinearAlignmentParameters
{
    LinearAlignmentParameters() { assertThatScoresAreValid(); }
    LinearAlignmentParameters(int32_t matchScore, int32_t mismatchScore, int32_t gapOpenScore, int32_t gapExtendScore)
        : matchScore(matchScore)
        , mismatchScore(mismatchScore)
        , gapOpenScore(gapOpenScore)
        , gapExtendScore(gapExtendScore)
    {
        assertThatScoresAreValid();
    }

    void assertThatScoresAreValid()
    {
        assert(0 <= matchScore && mismatchScore <= 0 && gapOpenScore <= 0 && gapExtendScore <= 0);
    }

    const int32_t matchScore = 5;
    const int32_t mismatchScore = -4;
    const int32_t gapOpenScore = -8;
    const int32_t gapExtendScore = -2;
};
