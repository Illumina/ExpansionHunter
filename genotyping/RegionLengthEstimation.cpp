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

#include "genotyping/RegionLengthEstimation.hh"

#include <boost/lexical_cast.hpp>
#include <boost/math/distributions.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

using boost::lexical_cast;
using boost::math::cdf;
using boost::math::poisson_distribution;
using std::round;
using std::vector;

namespace ehunter
{

// Uses the standard Lander-Waterman model to estimate length of a region. The confidence interval is computed using a
// generic parametric bootstrap procedure (note that a simpler implementation using Poisson mean CI is possible).
void estimateRegionLength(
    int readCount, int readLength, double depth, int& regionLength, int& lowerBound, int& upperBound)
{
    const double proportionOfReadsStartAtPosition = depth / readLength;
    // The length of sub-region where reads can start and still be fully within the region
    const int extensionLength = static_cast<int>(round(readCount / proportionOfReadsStartAtPosition));

    const int kSeed = 42;
    std::mt19937 numberGenerator(kSeed);

    // Model for the number of reads that fall within the region
    std::poisson_distribution<> poisson(readCount);

    vector<int> bootstrapSamples;
    const int kNumSamples = 10000;
    for (int sampleIndex = 0; sampleIndex < kNumSamples; ++sampleIndex)
    {
        const int sampledReadCount = poisson(numberGenerator);
        const int sampledExtensionLength = static_cast<int>(round(sampledReadCount / proportionOfReadsStartAtPosition));
        const int bootstrapSample = sampledExtensionLength - extensionLength;

        bootstrapSamples.push_back(bootstrapSample);
    }

    // Compute 2.5% and 97.5% quantiles
    std::sort(bootstrapSamples.begin(), bootstrapSamples.end());
    const int lowerQuantile = *(bootstrapSamples.begin() + static_cast<int>(bootstrapSamples.size() * 0.025));
    const int upperQuantile = *(bootstrapSamples.begin() + static_cast<int>(bootstrapSamples.size() * 0.975));

    regionLength = extensionLength + readLength;

    lowerBound = readLength;
    lowerBound += extensionLength - upperQuantile > 0 ? extensionLength - upperQuantile : 0;

    upperBound = readLength;
    upperBound += extensionLength - lowerQuantile > 0 ? extensionLength - lowerQuantile : 0;
}

}
