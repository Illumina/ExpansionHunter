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

namespace ehunter
{

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

}
