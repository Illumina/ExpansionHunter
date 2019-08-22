//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Konrad Scheffler <kscheffler@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "genotyping/ProbabilisticRepeatGenotyper.hh"

#include "gtest/gtest.h"

using namespace ehunter;

TEST(ProbabilisticallyGenotypingRepeats, ShortRepeatWithOneAllele_Genotyped)
{
    const int32_t repeatUnitLen = 1;
    const int32_t adjustedRegionSize = 100;
    const int32_t readLength = 3;
    const int32_t maxAlleleSize = 10;
    const double stutterPenalty = -2.0;
    const double randomBasePenalty = 0.0;
    const double mismapProb = 0.01;

    // Single read with single alignment spanning 2 repeat units
    ReadSummaryForStr summary(readLength);
    // alleleSize score readType clippedReadLength
    summary.addAlignment(StrAlignment(2, StrAlignment::Type::kSpanning, 5 * readLength, readLength));

    // Give genotyper two reads with identical summaries
    ProbabilisticRepeatGenotyper genotyper(
        AlleleCount::kOne, repeatUnitLen, adjustedRegionSize, readLength, maxAlleleSize, stutterPenalty,
        randomBasePenalty, mismapProb, { summary, summary });

    const optional<RepeatGenotype> expectedGenotype = RepeatGenotype(1, { 2 });

    const double credibleIntervalSize = 0.95;
    ASSERT_EQ(expectedGenotype, genotyper.genotypeRepeat(credibleIntervalSize));
}

TEST(ProbabilisticallyGenotypingRepeats, ShortRepeatWithTwoAlleles_Genotyped)
{
    const int32_t repeatUnitLen = 1;
    const int32_t adjustedRegionSize = 100;
    const int32_t readLength = 3;
    const int32_t maxAlleleSize = 10;
    const double stutterPenalty = -2.0;
    const double randomBasePenalty = 0.0;
    const double mismapProb = 0.01;

    // Two reads, each with one alignment spanning 2 and 3 repeat units respectively
    ReadSummaryForStr summary1(readLength);
    ReadSummaryForStr summary2(readLength);
    // alleleSize score readType clippedReadLength
    summary1.addAlignment(StrAlignment(2, StrAlignment::Type::kSpanning, 5 * readLength, readLength));
    summary2.addAlignment(StrAlignment(3, StrAlignment::Type::kSpanning, 5 * readLength, readLength));

    // Give genotyper two of each type of read
    ProbabilisticRepeatGenotyper genotyper(
        AlleleCount::kTwo, repeatUnitLen, adjustedRegionSize, readLength, maxAlleleSize, stutterPenalty,
        randomBasePenalty, mismapProb, { summary1, summary1, summary2, summary2 });

    RepeatGenotype expectedGenotype = RepeatGenotype(1, { 2, 3 });
    expectedGenotype.setShortAlleleSizeInUnitsCi(0, 3);
    expectedGenotype.setLongAlleleSizeInUnitsCi(0, 3);

    const double credibleIntervalSize = 0.95;
    ASSERT_EQ(expectedGenotype, *genotyper.genotypeRepeat(credibleIntervalSize));
}
