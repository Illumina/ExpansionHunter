//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
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

#include "genotyping/genotype.h"

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
