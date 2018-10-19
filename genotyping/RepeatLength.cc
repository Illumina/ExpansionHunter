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

#include "genotyping/RepeatLength.hh"

#include <boost/lexical_cast.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using boost::lexical_cast;
using boost::math::binomial_distribution;
using boost::math::cdf;
using std::vector;

// Given the observed IRR number, haplotype depth, and read length
// estimate repeat length (in nt) and the associated confidence
// interval.
void estimateRepeatLen(
    int32_t numIrrs, int32_t readLen, double hapDepth, int32_t& lenEstimate, int32_t& lowerBound, int32_t& upperBound)
{
    const double probReadStart = hapDepth / readLen;
    const int mlEstimate = static_cast<int>(std::round(numIrrs / probReadStart));

    const int32_t kSeed = 42;
    std::mt19937 gen(kSeed);

    // Perform mlEstimate trials with probability of succeed probReadStart.
    std::binomial_distribution<> binom(mlEstimate, probReadStart);

    vector<int32_t> bootstrapSamples;
    const int32_t kNumSamples = 10000;
    for (int32_t sampleIndex = 0; sampleIndex < kNumSamples; ++sampleIndex)
    {
        const int32_t sampledNumIrrs = binom(gen);
        const int32_t bootstrapSample = static_cast<int32_t>(std::round(sampledNumIrrs / probReadStart)) - mlEstimate;
        bootstrapSamples.push_back(bootstrapSample);
    }

    std::sort(bootstrapSamples.begin(), bootstrapSamples.end());

    // Compute 2.5% and 97.5% quantiles.
    const int32_t lowerQuantile = *(bootstrapSamples.begin() + static_cast<int32_t>(bootstrapSamples.size() * 0.025));
    const int32_t upperQuantile = *(bootstrapSamples.begin() + static_cast<int32_t>(bootstrapSamples.size() * 0.975));

    lenEstimate = mlEstimate + readLen;

    lowerBound = 0;
    if (mlEstimate - upperQuantile > 0)
    {
        lowerBound = static_cast<int32_t>(mlEstimate - upperQuantile);
    }
    lowerBound += readLen;

    assert(mlEstimate - lowerQuantile + readLen >= 0);
    upperBound = static_cast<int32_t>(mlEstimate - lowerQuantile + readLen);
}
