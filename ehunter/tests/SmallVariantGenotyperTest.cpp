//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Sai Chen <schen6@illumina.com>
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

#include "genotyping/SmallVariantGenotyper.hh"

#include <numeric>

#include "gtest/gtest.h"

#include "core/Common.hh"

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
