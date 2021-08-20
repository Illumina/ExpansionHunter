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

#include "genotyping/AlleleChecker.hh"

#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/log1p.hpp>

#include <numeric>

namespace ehunter
{

namespace
{
double poissonLogPmf(double lambda, double count)
{
    return count * log(lambda) - lambda - boost::math::lgamma(count + 1);
}

double logBeta(int a, int b) { return boost::math::lgamma(a) + boost::math::lgamma(b) - boost::math::lgamma(a + b); }

double logBinomCoef(int n, int k) { return -boost::math::log1p(n) - logBeta(n - k + 1, k + 1); }

double binomLogPmf(int n, double p, int count)
{
    return logBinomCoef(n, count) + count * log(p) + (n - count) * boost::math::log1p(-p);
}

}

AlleleCheckSummary AlleleChecker::check(double haplotypeDepth, int targetAlleleCount, int otherAlleleCount) const
{
    if (haplotypeDepth <= 0)
    {
        throw std::runtime_error("Haplotype depth must be positive");
    }

    if (targetAlleleCount < 0 || otherAlleleCount < 0)
    {
        throw std::runtime_error("Negative read counts are not allowed");
    }

    const int totalReadCount = targetAlleleCount + otherAlleleCount;
    const double ll0 = (totalReadCount > 0) ? binomLogPmf(totalReadCount, errorRate_, targetAlleleCount) : 0;
    const double ll1 = poissonLogPmf(haplotypeDepth, targetAlleleCount);

    AlleleStatus status = AlleleStatus::kUncertain;
    double logLikelihoodRatio = (ll1 - ll0) / log(10);
    if (logLikelihoodRatio < -log10(likelihoodRatioThreshold_))
    {
        status = AlleleStatus::kAbsent;
    }
    else if (logLikelihoodRatio > log10(likelihoodRatioThreshold_))
    {
        status = AlleleStatus::kPresent;
    }

    return AlleleCheckSummary(status, logLikelihoodRatio);
}

std::ostream& operator<<(std::ostream& out, AlleleStatus status)
{
    switch (status)
    {
    case AlleleStatus::kAbsent:
        out << "Absent";
        break;
    case AlleleStatus::kPresent:
        out << "Present";
        break;
    case AlleleStatus::kUncertain:
        out << "Uncertain";
        break;
    }

    return out;
}

}
