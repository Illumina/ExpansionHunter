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

double Poisson::Pmf(int32_t count) const {
  return boost::math::pdf(boost_poisson_, count);
}

double Poisson::ComputeProbabilityOfLessExtremeCounts(int32_t count) const {
  double total_probability = 0;

  for (int32_t current_count = count - 1; current_count >= 0; --current_count) {
    if (Pmf(current_count) > Pmf(count)) {
      total_probability += Pmf(current_count);
    } else {
      break;
    }
  }

  const int32_t kHardUpperLimit = 50000;
  for (int32_t current_count = count + 1; current_count != kHardUpperLimit;
       ++current_count) {
    if (Pmf(current_count) > Pmf(count)) {
      total_probability += Pmf(current_count);
    } else {
      break;
    }
  }

  return total_probability;
}

double Poisson::ComputeProbabilityOfCountsAsExtreme(int32_t count) const {
  return 1 - ComputeProbabilityOfLessExtremeCounts(count);
}

double ExpectedCountTest::Test(int32_t count) const {
  return poisson_.ComputeProbabilityOfCountsAsExtreme(count);
}