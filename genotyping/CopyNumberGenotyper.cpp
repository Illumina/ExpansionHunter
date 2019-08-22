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
#include "genotyping/CopyNumberGenotyper.hh"
#include <boost/math/distributions/normal.hpp>
#include <math.h>
#include <numeric>

using boost::optional;
using boost::math::normal_distribution;
using std::vector;

namespace ehunter
{

CopyNumberGenotyper::CopyNumberGenotyper(
    int maxCopyNumber, double depthScaleFactor, double standardDeviationOfCN2,
    const std::vector<double>& meanDepthValues, const std::vector<double>& priorCopyNumberFreq)
    : maxCopyNumber_(maxCopyNumber)
    , depthScaleFactor_(depthScaleFactor)
    , standardDeviationOfCN2_(standardDeviationOfCN2)
    , meanDepthValues_(meanDepthValues)
    , priorCopyNumberFreq_(priorCopyNumberFreq)
{
    if (maxCopyNumber_ + 1 != static_cast<int>(meanDepthValues_.size()))
    {
        throw std::runtime_error("Number of mean values is inconsistent with total copy number states.");
    }

    if (maxCopyNumber_ + 1 != static_cast<int>(priorCopyNumberFreq_.size()))
    {
        throw std::runtime_error("Number of prior frequencies is inconsistent with total copy number states.");
    }
}

boost::optional<int> CopyNumberGenotyper::genotype(double normalizedDepth) const
{
    const double adjustedDepth = normalizedDepth / depthScaleFactor_;
    std::vector<double> likelihoodOfAllCN;
    std::vector<double> pvalueOfAllCN;

    for (int currentGenotype = 0; currentGenotype != maxCopyNumber_ + 1; currentGenotype++)
    {
        std::pair<double, double> likelihoodAndPvalue = genotypeLikelihoodAndPvalue(currentGenotype, adjustedDepth);
        double currentLikelihood = likelihoodAndPvalue.first;
        double currentPvalue = likelihoodAndPvalue.second;
        likelihoodOfAllCN.emplace_back(currentLikelihood);
        pvalueOfAllCN.emplace_back(currentPvalue);
    }

    const std::pair<int, double> bestGenotypeAndPosterior = getBestGenotypeAndPosterior(likelihoodOfAllCN);
    const int bestGenotype = bestGenotypeAndPosterior.first;
    const double posteriorProbabilityOfBestGenotype = bestGenotypeAndPosterior.second;
    const bool posteriorCheck = posteriorProbabilityOfBestGenotype > posteriorProbabilityThreshold_;
    const bool pvalueCheck = pvalueOfAllCN[bestGenotype] > pvalueThreshold_;
    const optional<int> genotype = (posteriorCheck && pvalueCheck) ? bestGenotype : optional<int>();
    return genotype;
}

std::pair<int, double>
CopyNumberGenotyper::getBestGenotypeAndPosterior(const std::vector<double>& likelihoodOfAllCN) const
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

std::pair<double, double>
CopyNumberGenotyper::genotypeLikelihoodAndPvalue(int currentGenotype, double adjustedDepth) const
{
    assert(currentGenotype < (int)meanDepthValues_.size());
    assert(currentGenotype < (int)priorCopyNumberFreq_.size());

    const double meanValue = meanDepthValues_[currentGenotype];
    // standard deviation for each CN state can be computed from CN2 except for CN0 it is pre-defined
    const double standardDeviation
        = currentGenotype == 0 ? standardDeviationOfCN0_ : standardDeviationOfCN2_ * sqrt(double(currentGenotype) / 2);
    const double priorFreq = priorCopyNumberFreq_[currentGenotype];

    const normal_distribution<> cnDistribution(meanValue, standardDeviation);
    const double genotypeLikelihood = priorFreq * pdf(cnDistribution, adjustedDepth);
    const double cumulativeFrequency = cdf(cnDistribution, adjustedDepth);
    const double pvalue = std::min(cumulativeFrequency, 1 - cumulativeFrequency);

    std::pair<double, double> likelihoodAndPvalue;
    likelihoodAndPvalue.first = genotypeLikelihood;
    likelihoodAndPvalue.second = pvalue;

    return likelihoodAndPvalue;
}
}
