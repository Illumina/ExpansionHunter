//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "genotyping/RepeatGenotype.hh"

#include "gtest/gtest.h"

#include "common/common.h"

using std::vector;

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

TEST(SettingAlleleSizesCis, CiNotCoveringRepeatSize_ExpceptionThrown)
{
    RepeatGenotype genotype(3, { 2, 3 });

    EXPECT_ANY_THROW(genotype.setShortAlleleSizeInUnitsCi(0, 1));
    EXPECT_ANY_THROW(genotype.setLongAlleleSizeInUnitsCi(4, 5));
}

TEST(TestingHomozygosity, TypicalGenotypes_HomozygosityDetermined)
{
    EXPECT_TRUE(RepeatGenotype(3, { 2 }).isHomozygous());
    EXPECT_FALSE(RepeatGenotype(3, { 2, 3 }).isHomozygous());
    EXPECT_TRUE(RepeatGenotype(3, { 3, 3 }).isHomozygous());
}
