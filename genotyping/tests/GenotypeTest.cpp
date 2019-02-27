//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
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

#include "genotyping/genotype.h"

using namespace ehunter;

TEST(InitializingStrAlleles, TypicalAllele_Initialized)
{
    StrAllele allele(3, ReadType::kSpanning);
    EXPECT_EQ(3, allele.size());
    EXPECT_EQ(Interval(3, 3), allele.sizeRange());
    EXPECT_EQ(ReadType::kSpanning, allele.supportType());
}

TEST(InitializingStrAlleles, AlleleSupportedBySpanningReads_SizeRangeMustEqualToSize)
{
    StrAllele allele(3, ReadType::kSpanning);
    EXPECT_NO_THROW(allele.setSizeRange(3, 3));
    EXPECT_ANY_THROW(allele.setSizeRange(4, 4));
    EXPECT_ANY_THROW(allele.setSizeRange(2, 5));
}

TEST(InitializingStrAlleles, AlleleSupportedByFlankingOrRepeatReads_SizeRangeMustContainSize)
{
    StrAllele allele(10, ReadType::kFlanking);
    EXPECT_NO_THROW(allele.setSizeRange(5, 15));
    EXPECT_ANY_THROW(allele.setSizeRange(11, 12));
    EXPECT_ANY_THROW(allele.setSizeRange(8, 9));
}

// TEST(DefiningStrGenotypes, TypicalStrGenotype_Test)
//{
//    StrGenotype genotype = {{3, ReadType::kSpanning}, {50, ReadType::kRepeat}};

// EXPECT_EQ(3, genotype.shortAlleleSize());
// EXPECT_EQ(Interval(3, 3), genotype.shortAlleleCi());
// EXPECT_EQ(ReadType::kSpanning, genotype.shortAlleleSupportType());

// EXPECT_EQ(5, genotype.longAlleleSize());
// EXPECT_EQ(Interval(5, 5), genotype.longAlleleCi());
// EXPECT_EQ(ReadType::kRepeat, genotype.longAlleleSupportType());
//}
