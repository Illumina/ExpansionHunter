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

#include "include/repeat_length.h"

#include <boost/math/distributions.hpp>
using boost::math::cdf;
#include <boost/math/distributions/binomial.hpp>
using boost::math::binomial_distribution;
#include <boost/lexical_cast.hpp>
using boost::lexical_cast;

#include <cmath>
#include <iostream>
using std::cerr;
using std::endl;
#include <random>
#include <vector>
using std::vector;

// Given the observed IRR number, haplotype depth, and read length
// estimate repeat length (in nt) and the associated confidence
// interval.
void EstimateRepeatLen(const size_t num_irrs, const size_t read_len,
                       const double hap_depth, size_t& len_estimate,
                       size_t& lower_bound, size_t& upper_bound) {
  const double prob_read_start = hap_depth / read_len;
  const int ml_estimate =
      static_cast<int>(std::round(num_irrs / prob_read_start));

  std::random_device rd;
  std::mt19937 gen(rd());

  // Perform ml_estimate trials with probability of succeed prob_read_start.
  std::binomial_distribution<> binom(ml_estimate, prob_read_start);

  vector<int> bootstrap_samples;
  const size_t NUM_SAMPLES = 10000;
  for (size_t n = 0; n < NUM_SAMPLES; ++n) {
    const size_t sampled_num_irrs = binom(gen);
    const int bootstrap_sample =
        static_cast<int>(std::round(sampled_num_irrs / prob_read_start)) -
        ml_estimate;
    bootstrap_samples.push_back(bootstrap_sample);
  }

  std::sort(bootstrap_samples.begin(), bootstrap_samples.end());

  // Compute 2.5% and 97.5% quantiles.
  const int lower_quantile =
      *(bootstrap_samples.begin() +
        static_cast<int>(bootstrap_samples.size() * 0.025));
  const int upper_quantile =
      *(bootstrap_samples.begin() +
        static_cast<int>(bootstrap_samples.size() * 0.975));

  len_estimate = ml_estimate + read_len;

  lower_bound = 0;
  if (ml_estimate - upper_quantile > 0) {
    lower_bound = (size_t)(ml_estimate - upper_quantile);
  }
  lower_bound += read_len;

  assert(ml_estimate - lower_quantile + read_len >= 0);
  upper_bound = (size_t)(ml_estimate - lower_quantile + read_len);
}
