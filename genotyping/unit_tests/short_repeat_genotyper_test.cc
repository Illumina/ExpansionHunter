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

#include "genotyping/short_repeat_genotyper.h"

#include <array>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "common/common.h"

using std::vector;
using std::string;
using std::map;
using std::cerr;
using std::endl;
using std::array;

using namespace ehunter;

TEST(CalculateMoleculeProportions, TypicalHaplotypeProportionsCalculated) {
  const int num_units_haplotype = 2;
  const int max_num_units_in_read = 25;
  const double prop_correct_molecules = 0.97;
  Allele hap(num_units_haplotype, max_num_units_in_read,
             prop_correct_molecules);

  EXPECT_DOUBLE_EQ(2.2885056508333023e-08, hap.propMolecules(25));
  EXPECT_DOUBLE_EQ(0.97087262363952287, hap.propMoleculesShorterThan(3));
  EXPECT_DOUBLE_EQ(0.029127376360477131, hap.propMoleculesAtLeast(3));
}

TEST(CalcFlankingLoglik, TypicalFlankingReadsLoglikelihoodsCalculated) {
  ShortRepeatGenotyper genotype(25, 0.97, 20.0, 150, 2, 3);
  EXPECT_DOUBLE_EQ(-2.0300033341853156, genotype.CalcFlankingLoglik(2));
  EXPECT_DOUBLE_EQ(-19.607697373350305, genotype.CalcFlankingLoglik(25));
}

TEST(CalcSpanningLoglik, TypicalSpanningReadsLoglikelihoodsCalcualted) {
  ShortRepeatGenotyper genotype(25, 0.97, 20.0, 150, 2, 3);
  EXPECT_DOUBLE_EQ(-2.7385082705573418, genotype.CalcSpanningLoglik(3));
  EXPECT_DOUBLE_EQ(-6.2450661678773223, genotype.CalcSpanningLoglik(4));
}

TEST(CalcGenotypeLoglik, ShortGenotypesLoglikelihoodsCalculated) {
  const map<int, int> flanking_size_counts = {{1, 2}, {2, 3}, {10, 1}};
  const map<int, int> spanning_size_counts = {{3, 4}, {5, 5}};
  vector<AlleleSupport> genotype_support;

  ShortRepeatGenotyper genotype_3_5(25, 0.97, 25.0, 150, 3, 5);
  EXPECT_DOUBLE_EQ(
      -48.468337669679954,
      genotype_3_5.CalcLogLik(flanking_size_counts, spanning_size_counts,
                              genotype_support));

  const vector<AlleleSupport> expected_3_5_support = {{4, 5, 0}, {5, 5, 0}};
  EXPECT_EQ(expected_3_5_support, genotype_support);

  ShortRepeatGenotyper genotype_3_10(25, 0.97, 25.0, 150, 3, 10);
  EXPECT_DOUBLE_EQ(
      -69.444360064064853,
      genotype_3_10.CalcLogLik(flanking_size_counts, spanning_size_counts,
                               genotype_support));
  const vector<AlleleSupport> expected_3_10_support = {{4, 5, 0}, {0, 6, 0}};
  EXPECT_EQ(expected_3_10_support, genotype_support);

  ShortRepeatGenotyper genotype_10_10(25, 0.97, 25.0, 150, 10, 10);
  EXPECT_DOUBLE_EQ(
      -185.24122167420646,
      genotype_10_10.CalcLogLik(flanking_size_counts, spanning_size_counts,
                                genotype_support));
  const vector<AlleleSupport> expected_10_10_support = {{0, 6, 0}, {0, 6, 0}};
  EXPECT_EQ(expected_10_10_support, genotype_support);
}

TEST(CalcGenotypeLoglik, LongGenotypesLoglikelihoodsCalculated) {
  const map<int, int> flanking_size_counts = {{1, 2}, {2, 3}, {10, 1}};
  const map<int, int> spanning_size_counts = {{3, 4}, {5, 5}};
  vector<AlleleSupport> genotype_support;

  ShortRepeatGenotyper genotype_3_5(25, 0.97, 25.0, 150, 3, 5);
  EXPECT_DOUBLE_EQ(
      -48.468337669679954,
      genotype_3_5.CalcLogLik(flanking_size_counts, spanning_size_counts,
                              genotype_support));
}

TEST(CalcDiploidGenotypeLoglik, TypicalGenotypeLoglikelihoodsCalculated) {
  const map<int, int> flanking_size_counts = {{1, 2}, {2, 3}, {25, 10}};
  const map<int, int> spanning_size_counts = {{5, 5}};
  vector<AlleleSupport> genotype_support;

  ShortRepeatGenotyper diploid_genotype(25, 0.97, 25.0, 150, 5, 25);
  EXPECT_DOUBLE_EQ(
      -34.260255045398637,
      diploid_genotype.CalcLogLik(flanking_size_counts, spanning_size_counts,
                                  genotype_support));

  const vector<AlleleSupport> expected_5_25_support = {{5, 5, 0}, {0, 5, 10}};
  EXPECT_EQ(expected_5_25_support, genotype_support);
}

TEST(GenotypeStr, TypicalDiploidStrReturnsGenotype) {
  const map<int, int> flanking_size_counts = {{1, 2}, {2, 3}, {10, 1}};
  const map<int, int> spanning_size_counts = {{3, 4}, {5, 5}};

  const int max_num_units_in_read = 25;
  const double prop_correct_molecules = 0.97;
  const double hap_depth = 25.0;
  const int read_len = 150;

  vector<RepeatAllele> haplotype_candidates;
  for (int i = 0; i != 26; ++i) {
    haplotype_candidates.push_back(RepeatAllele(i, -1, ReadType::kSpanning));
  }

  RepeatGenotype genotype;

  GenotypeShortRepeat(max_num_units_in_read, prop_correct_molecules, hap_depth,
                      read_len, haplotype_candidates, flanking_size_counts,
                      spanning_size_counts, GenotypeType::kDiploid, genotype);

  const vector<RepeatAllele> expected_genotype = {
      RepeatAllele(3, ReadType::kSpanning, AlleleSupport(4, 5, 0)),
      RepeatAllele(5, ReadType::kSpanning, AlleleSupport(5, 5, 0))};
  EXPECT_EQ(expected_genotype, genotype);
}

TEST(GenotypeStr, TypicalHaploidStrReturnsGenotype) {
  const map<int, int> flanking_size_counts = {{1, 2}, {2, 3}, {10, 1}};
  const map<int, int> spanning_size_counts = {{3, 4}, {5, 5}};

  const int max_num_units_in_read = 25;
  const double prop_correct_molecules = 0.97;
  const double hap_depth = 25.0;
  const int read_len = 150;

  vector<RepeatAllele> haplotype_candidates;
  for (int i = 0; i != 26; ++i) {
    haplotype_candidates.push_back(RepeatAllele(i, -1, ReadType::kSpanning));
  }

  RepeatGenotype genotype;
  GenotypeShortRepeat(max_num_units_in_read, prop_correct_molecules, hap_depth,
                      read_len, haplotype_candidates, flanking_size_counts,
                      spanning_size_counts, GenotypeType::kHaploid, genotype);

  const vector<RepeatAllele> expected_genotype = {
      RepeatAllele(5, ReadType::kSpanning, AlleleSupport(5, 5, 0))};

  EXPECT_EQ(expected_genotype, genotype);
}

TEST(GenotypeStr, ExpandedHaploidStrGenotyped) {
  const map<int, int> flanking_size_counts = {{1, 2}, {2, 3}, {10, 1}, {25, 8}};
  const map<int, int> spanning_size_counts = {{3, 1}, {5, 1}};

  const int max_num_units_in_read = 25;
  const double prop_correct_molecules = 0.97;
  const double hap_depth = 25.0;
  const int read_len = 150;

  vector<RepeatAllele> haplotype_candidates;
  for (int i = 0; i != 26; ++i) {
    haplotype_candidates.push_back(RepeatAllele(i, -1, ReadType::kSpanning));
  }

  RepeatGenotype genotype;
  GenotypeShortRepeat(max_num_units_in_read, prop_correct_molecules, hap_depth,
                      read_len, haplotype_candidates, flanking_size_counts,
                      spanning_size_counts, GenotypeType::kHaploid, genotype);

  const vector<RepeatAllele> expected_genotype = {
      RepeatAllele(25, ReadType::kSpanning, AlleleSupport(0, 6, 8))};

  EXPECT_EQ(expected_genotype, genotype);
}
