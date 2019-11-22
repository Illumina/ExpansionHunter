//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Xiao Chen <xchen2@illumina.com>
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
//

#include "genotyping/SmallVariantCopyNumberGenotyper.hh"
#include "gtest/gtest.h"
#include <numeric>
#include <vector>

#include "common/Common.hh"

using boost::optional;
using namespace ehunter;

TEST(CopyNumberGenotyping, RegularGenotype)
{
    {
        int totalCopyNumber = 4;
        SmallVariantCopyNumberGenotyper genotyper(totalCopyNumber);

        auto bestGenotype = *genotyper.genotype(0, 40, 3);
        EXPECT_EQ(0, bestGenotype.first);

        bestGenotype = *genotyper.genotype(10, 30, 3);
        EXPECT_EQ(1, bestGenotype.first);

        bestGenotype = *genotyper.genotype(20, 20, 3);
        EXPECT_EQ(2, bestGenotype.first);

        bestGenotype = *genotyper.genotype(30, 10, 3);
        EXPECT_EQ(3, bestGenotype.first);

        bestGenotype = *genotyper.genotype(39, 1, 3);
        EXPECT_EQ(4, bestGenotype.first);
    }
    
    {
        int totalCopyNumber = 3;
        SmallVariantCopyNumberGenotyper genotyper(totalCopyNumber);
        auto bestGenotype = *genotyper.genotype(30, 30, 3);
        EXPECT_EQ(2, bestGenotype.first);
        EXPECT_NEAR(0.689, bestGenotype.second, 1e-3);
    }

    {
        int totalCopyNumber = 3;
        SmallVariantCopyNumberGenotyper genotyper(totalCopyNumber);
        auto bestGenotype = *genotyper.genotype(3, 10, 4);
        EXPECT_EQ(0, bestGenotype.first);
        bestGenotype = *genotyper.genotype(3, 10, 2);
        EXPECT_EQ(1, bestGenotype.first);
    }
    
}
