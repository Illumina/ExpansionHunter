//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "genotyping/RepeatGenotyper.hh"

#include "gtest/gtest.h"

#include "common/Common.hh"

using std::map;
using std::vector;

using namespace ehunter;

TEST(CountingInrepeatReads, HaploidExpansion_IrrsCounted)
{
    CountTable countsOfFlankingReads(map<int, int>(
        { { 1, 3 }, { 2, 3 }, { 7, 1 }, { 11, 1 }, { 18, 1 }, { 20, 1 }, { 21, 1 }, { 33, 1 }, { 44, 1 } }));
    CountTable countsOfInrepeatReads(map<int, int>({ { 43, 1 }, { 45, 6 }, { 46, 1 }, { 47, 2 }, { 48, 1 } }));

    const int maxNumUnitsInRead = 50;
    EXPECT_EQ(10, countFullLengthRepeatReads(maxNumUnitsInRead, countsOfFlankingReads, countsOfInrepeatReads));
}

TEST(CountingInrepeatReads, HaploidNormal_IrrsCounted)
{
    CountTable countsOfFlankingReads(map<int, int>(
        { { 1, 3 }, { 2, 3 }, { 7, 1 }, { 11, 1 }, { 18, 1 }, { 20, 1 }, { 21, 1 }, { 33, 1 }, { 44, 1 } }));
    CountTable countsOfInrepeatReads(map<int, int>({ { 46, 1 }, { 47, 1 }, { 48, 1 } }));

    const int maxNumUnitsInRead = 50;
    EXPECT_EQ(1, countFullLengthRepeatReads(maxNumUnitsInRead, countsOfFlankingReads, countsOfInrepeatReads));
}

TEST(GenotypeExtension, NoFlankingReadsWhenSomeAreExpected_ExtensionAborted)
{
    // haplotypeDepth, expectedAlleleCount, repeatUnitLen, maxNumUnitsInRead, propCorrectMolecules, ...
    RepeatGenotyper genotyper(20, AlleleCount::kTwo, 10, 15, 0.8, CountTable(), CountTable(), CountTable(), 0);
    RepeatGenotype genotype(10, { 2, 10 });
    RepeatGenotype expectedGenotype(genotype);

    genotyper.extendGenotypeWhenBothAllelesAreFlanking(genotype);
    EXPECT_EQ(expectedGenotype, genotype);

    genotyper.extendGenotypeWhenOneAlleleIsFlanking(genotype);
    EXPECT_EQ(expectedGenotype, genotype);
}
