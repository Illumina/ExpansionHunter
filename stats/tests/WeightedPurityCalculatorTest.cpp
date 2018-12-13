//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "stats/WeightedPurityCalculator.hh"

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
