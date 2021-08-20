//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
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
//

#include "core/WeightedPurityCalculator.hh"

#include <vector>

#include "gmock/gmock.h"

using std::string;
using std::vector;
using testing::DoubleNear;

using namespace ehunter;

TEST(CalculatingWeightedPurityScore, PerfectRepeat_Calculated)
{
    const string sequence = "GGCCCCGGCCCC";
    WeightedPurityCalculator wpCalculator("GGCCGG");
    EXPECT_THAT(wpCalculator.score(sequence), DoubleNear(1.0, 0.005));
}

TEST(CalculatingWeightedPurityScore, ImperfectRepeat_Calculated)
{
    WeightedPurityCalculator wpCalculator("AACCCC");
    EXPECT_THAT(wpCalculator.score("ACCCCAACCCCAACCCCAACCCCAACCCCAACCCCA"), DoubleNear(1.0, 0.005));
    EXPECT_THAT(wpCalculator.score("tCCCCttCCCCttCCCCttCCCCtTCCCCttCCCCT"), DoubleNear(0.75, 0.005));
}
