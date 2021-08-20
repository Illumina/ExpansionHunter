//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
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

#include <iostream>

#include "core/Common.hh"

namespace ehunter
{

enum class AlleleStatus
{
    kPresent,
    kAbsent,
    kUncertain
};

/// Results from the AlleleChecker on one allele
struct AlleleCheckSummary
{
    AlleleCheckSummary(AlleleStatus status, double logLikelihoodRatio)
        : status(status)
        , logLikelihoodRatio(logLikelihoodRatio)
    {
    }
    AlleleStatus status;
    double logLikelihoodRatio; // Log10(LR) for the allele being present
};

/*
 * Genotyper checking for presence (>= 1 allele) of a given 'Key' allele
 */
class AlleleChecker
{
public:
    AlleleChecker(double errorRate, double llrThreshold)
        : errorRate_(errorRate)
        , likelihoodRatioThreshold_(llrThreshold)
    {
        if (errorRate <= 0 || errorRate >= 1)
        {
            throw std::runtime_error("Error rate must be positive and less than 1");
        }
        if (llrThreshold < 0)
        {
            throw std::runtime_error("Likelihood Ratio threshold must be positive");
        }
    }

    AlleleCheckSummary check(double haplotypeDepth, int targetAlleleCount, int otherAlleleCount) const;

private:
    // Rate of 'false' key-allele observations
    double errorRate_;
    // If the likelihood ratio threshold in favor of presence or absence
    // is not at least this strong, return no call.
    double likelihoodRatioThreshold_;
};

std::ostream& operator<<(std::ostream& out, AlleleStatus status);

}
