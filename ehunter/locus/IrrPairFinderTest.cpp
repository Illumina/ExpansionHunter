//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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

#include "locus/IrrPairFinder.hh"

#include "gmock/gmock.h"

using namespace ehunter;
using namespace locus;

TEST(CreatingIrrPairFinder, TypicalArguments_Created)
{
    IrrPairFinder finder("CGG");
    ASSERT_EQ("CGG", finder.targetMotif());
}

TEST(CheckingForIrrPairs, PairOfPerfectIrrs_CheckPassed)
{
    IrrPairFinder finder("CGG");
    ASSERT_TRUE(finder.check("CGGCGGCG", "GGCGGC"));
}

TEST(CheckingForIrrPairs, PairOfIrrsWithWrongMotif_CheckFailed)
{
    IrrPairFinder finder("ATA");
    ASSERT_FALSE(finder.check("CGGCGGCG", "GGCGGC"));
}

TEST(CheckingForIrrPairs, NotIrrPair_CheckFailed)
{
    IrrPairFinder finder("CGG");
    ASSERT_FALSE(finder.check("ATACT", "GGCGGC"));
    ASSERT_FALSE(finder.check("GGCGGC", "ATACT"));
    ASSERT_FALSE(finder.check("ATACT", "ATACT"));
}
