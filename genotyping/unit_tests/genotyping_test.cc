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

#include "genotyping/genotyping.h"

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "gtest/gtest.h"

using std::vector;
using std::string;
using std::map;
using std::cerr;
using std::endl;

TEST(CalculateMoleculeProportions, TypicalHaplotypeProportionsCalculated)
{
  const int num_units_haplotype = 2;
  const int max_num_units_in_read = 25;
  const double prop_correct_molecules = 0.97;
  StrHaplotype hap(num_units_haplotype, max_num_units_in_read, prop_correct_molecules);

  EXPECT_DOUBLE_EQ(2.2885056508333023e-08, hap.propMolecules(25));
  EXPECT_DOUBLE_EQ(0.97087262363952287, hap.propMoleculesShorterThan(3));
  EXPECT_DOUBLE_EQ(0.029127376360477131, hap.propMoleculesAtLeast(3));
}

TEST(CalcFlankingLoglik, TypicalFlankingReadsLoglikelihoodsCalculated)
{
  StrGenotype genotype(25, 0.97, 20.0, 150, 2, 3);
  EXPECT_DOUBLE_EQ(-2.0300033341853156, genotype.calcFlankingLoglik(2));
  EXPECT_DOUBLE_EQ(-19.607697373350305, genotype.calcFlankingLoglik(25));
}

TEST(CalcSpanningLoglik, TypicalSpanningReadsLoglikelihoodsCalcualted)
{
  StrGenotype genotype(25, 0.97, 20.0, 150, 2, 3);
  EXPECT_DOUBLE_EQ(-2.7385082705573418, genotype.calcSpanningLoglik(3));
  EXPECT_DOUBLE_EQ(-6.2450661678773223, genotype.calcSpanningLoglik(4));
}

TEST(CalcGenotypeLoglik, TypicalGenotypesLoglikelihoodsCalculated)
{
  const map<int, int> flanking_size_counts = {{1, 2}, {2, 3}, {10, 1}};
  const map<int, int> spanning_size_counts = {{3, 4}, {5, 5}};

  StrGenotype genotype_3_5(25, 0.97, 25.0, 150, 3, 5);
  EXPECT_DOUBLE_EQ(-48.468337669679954, genotype_3_5.calcLogLik(flanking_size_counts, spanning_size_counts));

  StrGenotype genotype_3_10(25, 0.97, 25.0, 150, 3, 10);
  EXPECT_DOUBLE_EQ(-69.444360064064853, genotype_3_10.calcLogLik(flanking_size_counts, spanning_size_counts));

  StrGenotype genotype_10_10(25, 0.97, 25.0, 150, 10, 10);
  EXPECT_DOUBLE_EQ(-185.24122167420646, genotype_10_10.calcLogLik(flanking_size_counts, spanning_size_counts));
}