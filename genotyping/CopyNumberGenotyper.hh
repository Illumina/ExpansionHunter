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

#pragma once

#include <boost/optional.hpp>
#include <vector>

#include "common/Common.hh"

namespace ehunter
{

class CopyNumberGenotyper
{
public:
    CopyNumberGenotyper(
        int maxCopyNumber, double depthScaleFactor, double standardDeviationOfCN2,
        const std::vector<double>& meanDepthValues, const std::vector<double>& priorCopyNumberFreq);

    boost::optional<int> genotype(double normalizedDepth) const;

    // return genotype likelihood and p-value of the given genotype
    std::pair<double, double> genotypeLikelihoodAndPvalue(int currentGenotype, double adjustedDepth) const;
    /**
     * return the best genotype and its posterior probability
     * @param likelihoodOfAllCN vector containing likelihoods of each copy number state
     */
    std::pair<int, double> getBestGenotypeAndPosterior(const std::vector<double>& likelihoodOfAllCN) const;

private:
    int maxCopyNumber_;
    double depthScaleFactor_;
    // standard deviation for copy number 2
    double standardDeviationOfCN2_;
    // vector containing mean values for all possible copy number states
    std::vector<double> meanDepthValues_;
    // vector containing prior frequencies for all possible copy number states
    std::vector<double> priorCopyNumberFreq_;
    // hard-coded parameters
    // standard deviation for copy number 0
    double standardDeviationOfCN0_ = 0.032;
    // posterior probability threshold for calling copy number genotype
    double posteriorProbabilityThreshold_ = 0.95;
    // p-value threshold for calling copy number genotype
    double pvalueThreshold_ = 1e-3;
};
}
