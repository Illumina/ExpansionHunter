//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Xiao Chen <xchen2@illumina.com>
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "genotyping/SmallVariantCopyNumberGenotyper.hh"

#include <boost/math/distributions/poisson.hpp>
#include <numeric>
#include <iostream>

using boost::optional;
using boost::math::poisson_distribution;
using std::vector;

namespace ehunter
{

static std::pair<int, double> getBestGenotypeAndPosterior(const std::vector<double>& likelihoodOfAllCN)
{
    double sumOfLikelihood = 0;
    for (double likelihoodValue : likelihoodOfAllCN)
    {
        sumOfLikelihood += likelihoodValue;
    }

    assert(!likelihoodOfAllCN.empty());
    auto maxElement = max_element(likelihoodOfAllCN.begin(), likelihoodOfAllCN.end());
    const double maxLikelihood = *maxElement;
    std::pair<int, double> bestGenotypeAndPosterior;
    bestGenotypeAndPosterior.first = distance(likelihoodOfAllCN.begin(), maxElement);
    bestGenotypeAndPosterior.second = maxLikelihood / sumOfLikelihood;

    return bestGenotypeAndPosterior;
}

double SmallVariantCopyNumberGenotyper::genotypeLikelihood(
    int totalCopyNumber, int currentCopyNumber, int variantCount, int nonvariantCount) const
{
    int totalCount = variantCount + nonvariantCount;
    double depthExpected;
    if (currentCopyNumber == 0)
    {
        depthExpected = (errorRate_ / 3) * totalCount;
    }
    else if (currentCopyNumber == totalCopyNumber)
    {
        depthExpected = totalCount - errorRate_ * totalCount;
    }
    else
    {
        depthExpected = totalCount * currentCopyNumber / totalCopyNumber;
    }

    const poisson_distribution<> countDistribution(depthExpected);
    double likelihood = (variantCount <= nonvariantCount) ? pdf(countDistribution, variantCount)
                                                          : pdf(countDistribution, nonvariantCount);
    return likelihood;
}

boost::optional<std::pair<int, double>>
SmallVariantCopyNumberGenotyper::genotype(int variantCount, int nonvariantCount, int minReadSupport) const
{

    if (variantCount < 0 || nonvariantCount < 0)
    {
        throw std::runtime_error(
            "Invalid read counts: " + std::to_string(variantCount) + " " + std::to_string(nonvariantCount));
    }

    const int totalReadCount = variantCount + nonvariantCount;
    if (totalReadCount == 0) // missing genotype
    {
        return optional<std::pair<int, double>>();
    }

    std::vector<double> likelihoodOfAllCN;

    for (int currentCopyNumber = 0; currentCopyNumber != totalCopyNumber_ + 1; currentCopyNumber++)
    {
        const double currentLikelihood
            = genotypeLikelihood(totalCopyNumber_, currentCopyNumber, variantCount, nonvariantCount);
        likelihoodOfAllCN.emplace_back(currentLikelihood);
    }

    if (variantCount > nonvariantCount)
    {
        std::reverse(likelihoodOfAllCN.begin(), likelihoodOfAllCN.end());
    }

    const std::pair<int, double> bestGenotypeAndPosterior = getBestGenotypeAndPosterior(likelihoodOfAllCN);

    if (bestGenotypeAndPosterior.first != 0 && variantCount <= minReadSupport && nonvariantCount >= minReadSupport)
    {
        return std::pair<int, double>(0, 1);
    }

    return bestGenotypeAndPosterior;
}
}
