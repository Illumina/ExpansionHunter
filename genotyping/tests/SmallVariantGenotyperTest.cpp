//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Sai Chen <schen6@illumina.com>
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

#include "genotyping/SmallVariantGenotyper.hh"

#include <numeric>

#include "gtest/gtest.h"

#include "common/Common.hh"

using namespace ehunter;

TEST(SmallVariantGenotyper, ThrowsWithIllegalParameter)
{
    SmallVariantGenotyper triploidGenotyper(30.0, (AlleleCount)3);
    ASSERT_ANY_THROW(triploidGenotyper.genotype(20, 20));

    SmallVariantGenotyper negCountGenotyper(30.0, (AlleleCount)2);
    ASSERT_ANY_THROW(negCountGenotyper.genotype(-1, 20));
}

TEST(SmallVariantGenotyper, RegularGenotye)
{
    SmallVariantGenotyper genotyper(30.0, (AlleleCount)2);

    SmallVariantGenotype ref_genotype(AlleleType::kRef, AlleleType::kRef);
    SmallVariantGenotype gt0 = *genotyper.genotype(20, 1);
    EXPECT_EQ(ref_genotype, gt0);

    SmallVariantGenotype het_genotype(AlleleType::kRef, AlleleType::kAlt);
    SmallVariantGenotype gt1 = *genotyper.genotype(20, 19);
    EXPECT_EQ(het_genotype, gt1);

    SmallVariantGenotype alt_genotype(AlleleType::kAlt, AlleleType::kAlt);
    SmallVariantGenotype gt2 = *genotyper.genotype(1, 20);
    EXPECT_EQ(alt_genotype, gt2);
}