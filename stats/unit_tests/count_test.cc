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

#include "stats/counts.h"

#include "gmock/gmock.h"

#include <iostream>
#include <string>
#include <vector>

using std::cerr;
using std::endl;
using std::string;
using std::vector;

TEST(Poisson, CalculatesProbabilitiesOfTypicalObservations) {
  const double rate = 2.5;
  Poisson poisson(rate);

  const double max_err = 0.01;
  EXPECT_NEAR(0.082, poisson.Pmf(0), max_err);
  EXPECT_NEAR(0.205, poisson.Pmf(1), max_err);
  EXPECT_NEAR(0.257, poisson.Pmf(2), max_err);
  EXPECT_NEAR(0.213, poisson.Pmf(3), max_err);
  EXPECT_NEAR(0.133, poisson.Pmf(4), max_err);
  EXPECT_NEAR(0.067, poisson.Pmf(5), max_err);
  EXPECT_NEAR(0.028, poisson.Pmf(6), max_err);
  EXPECT_NEAR(0.010, poisson.Pmf(7), max_err);
}

TEST(Poisson, ComputesTotalProbabilityOfCountsLessExtremeThanGiven) {
  const double rate = 3.5;
  Poisson poisson(rate);

  const double max_err = 0.01;
  EXPECT_NEAR(0.943063694486, poisson.ComputeProbabilityOfLessExtremeCounts(0),
              max_err);
  EXPECT_NEAR(0.721725327695, poisson.ComputeProbabilityOfLessExtremeCounts(1),
              max_err);
  EXPECT_NEAR(0.404597754447, poisson.ComputeProbabilityOfLessExtremeCounts(2),
              max_err);
  EXPECT_NEAR(0.000000000000, poisson.ComputeProbabilityOfLessExtremeCounts(3),
              max_err);
  EXPECT_NEAR(0.215785469039, poisson.ComputeProbabilityOfLessExtremeCounts(4),
              max_err);
  EXPECT_NEAR(0.589556727909, poisson.ComputeProbabilityOfLessExtremeCounts(5),
              max_err);
  EXPECT_NEAR(0.827416169673, poisson.ComputeProbabilityOfLessExtremeCounts(6),
              max_err);
  EXPECT_NEAR(0.904514519549, poisson.ComputeProbabilityOfLessExtremeCounts(7),
              max_err);
  EXPECT_NEAR(0.973261077909, poisson.ComputeProbabilityOfLessExtremeCounts(8),
              max_err);
  EXPECT_NEAR(0.990126341944, poisson.ComputeProbabilityOfLessExtremeCounts(9),
              max_err);
  EXPECT_NEAR(0.996685055735, poisson.ComputeProbabilityOfLessExtremeCounts(10),
              max_err);
}

TEST(Poisson, ComputeProbabilityOfCountsAsExtremeAsGiven) {
  const double rate = 3.5;
  Poisson poisson(rate);

  const double max_err = 0.01;
  EXPECT_NEAR(1 - 0.943063694486,
              poisson.ComputeProbabilityOfCountsAsExtreme(0), max_err);
  EXPECT_NEAR(1 - 0.589556727909,
              poisson.ComputeProbabilityOfCountsAsExtreme(5), max_err);
  EXPECT_NEAR(1 - 0.996685055735,
              poisson.ComputeProbabilityOfCountsAsExtreme(10), max_err);
}

TEST(ExpectedCountTest, CalculatesPvaluesForTypicalCounts) {
  const int32_t expected_count = 30;
  ExpectedCountTest count_test(expected_count);
  const double max_err = 0.00001;
  EXPECT_NEAR(1.6842083283563625e-13, count_test.Test(0), max_err);
  EXPECT_NEAR(4.947383377140735e-05, count_test.Test(10), max_err);
  EXPECT_NEAR(0.067594242342400968, count_test.Test(20), max_err);
  EXPECT_NEAR(1.000000000000000000, count_test.Test(30), max_err);
  EXPECT_NEAR(0.081538139623059447, count_test.Test(40), max_err);
  EXPECT_NEAR(0.00068661378232848236, count_test.Test(50), max_err);
  EXPECT_NEAR(1.4486113849576654e-06, count_test.Test(60), max_err);
}
