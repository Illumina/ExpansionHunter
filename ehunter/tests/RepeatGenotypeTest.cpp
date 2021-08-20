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

#include "genotyping/RepeatGenotype.hh"

#include "gtest/gtest.h"

#include "core/Common.hh"

using std::vector;

using namespace ehunter;

TEST(InitializingRepeatGenotype, HaploidGenotype_Initialized)
{
    const int32_t repeatUnitLen = 3;
    RepeatGenotype genotype(repeatUnitLen, { 2 });

    EXPECT_EQ(1, genotype.numAlleles());
    EXPECT_EQ(2, genotype.shortAlleleSizeInUnits());
    EXPECT_EQ(2, genotype.longAlleleSizeInUnits());
}

TEST(InitializingRepeatGenotype, DiploidGenotype_Initialized)
{
    RepeatGenotype genotype(3, { 2, 3 });

    EXPECT_EQ(2, genotype.numAlleles());
    EXPECT_EQ(2, genotype.shortAlleleSizeInUnits());
    EXPECT_EQ(3, genotype.longAlleleSizeInUnits());
}

TEST(ExtractingAlleleSizesInBases, DiploidGenotype_SizesExtracted)
{
    RepeatGenotype genotype(3, { 2, 3 });

    EXPECT_EQ(6, genotype.shortAlleleSizeInBp());
    EXPECT_EQ(9, genotype.longAlleleSizeInBp());
}

TEST(InitializingRepeatGenotype, NetherDiploidNorHaploidGenotype_ExceptionThrown)
{
    EXPECT_ANY_THROW(RepeatGenotype genotype(3, {}));
    EXPECT_ANY_THROW(RepeatGenotype genotype(3, { 1, 2, 3 }));
}

TEST(InitializingRepeatGenotype, UnorderedAlleleSizes_ExceptionThrown)
{
    EXPECT_ANY_THROW(RepeatGenotype genotype(3, { 5, 2 }));
}

TEST(SettingAlleleSizesCis, TypicalGenotype_CiSet)
{
    RepeatGenotype genotype(3, { 2, 3 });

    genotype.setShortAlleleSizeInUnitsCi(1, 5);
    genotype.setLongAlleleSizeInUnitsCi(2, 8);

    EXPECT_EQ(NumericInterval(1, 5), genotype.shortAlleleSizeInUnitsCi());
    EXPECT_EQ(NumericInterval(2, 8), genotype.longAlleleSizeInUnitsCi());
}

TEST(SettingAlleleSizesCis, CiNotCoveringRepeatSize_CiSizeExtended)
{
    RepeatGenotype genotype(3, { 2, 3 });
    genotype.setShortAlleleSizeInUnitsCi(0, 1);
    genotype.setLongAlleleSizeInUnitsCi(4, 5);

    EXPECT_EQ(NumericInterval(0, 2), genotype.shortAlleleSizeInUnitsCi());
    EXPECT_EQ(NumericInterval(3, 5), genotype.longAlleleSizeInUnitsCi());
}

TEST(TestingHomozygosity, TypicalGenotypes_HomozygosityDetermined)
{
    EXPECT_TRUE(RepeatGenotype(3, { 2 }).isHomozygous());
    EXPECT_FALSE(RepeatGenotype(3, { 2, 3 }).isHomozygous());
    EXPECT_TRUE(RepeatGenotype(3, { 3, 3 }).isHomozygous());
}
