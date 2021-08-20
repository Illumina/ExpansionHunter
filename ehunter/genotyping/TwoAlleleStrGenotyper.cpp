//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Konrad Scheffler <kscheffler@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

#include "genotyping/TwoAlleleStrGenotyper.hh"

#include <stack>

#include "core/LogSum.hh"

using std::unordered_set;

namespace ehunter
{
namespace strgt
{

struct Ci
{
    Ci(int begin, int end)
        : begin(begin)
        , end(end)
    {
    }
    int begin;
    int end;
};

struct CiAndLoglik
{
    CiAndLoglik(int startSize, int endSize, double loglik)
        : startSize(startSize)
        , endSize(endSize)
        , loglik(loglik)
    {
    }
    int startSize;
    int endSize;
    double loglik;
};

using LikelihoodFn = double (TwoAlleleGenotyper::*)(int, int);
static Ci getCiAlongX(int bestX, int bestY, TwoAlleleGenotyper* genotyper, LikelihoodFn likelihood)
{
    int xFrom = bestX;
    int xTo = bestX;
    int yFrom = bestY;
    int yTo = bestY;

    std::stack<CiAndLoglik> ciCandidates;
    double topGtLoglik = (genotyper->*likelihood)(bestX, bestY);
    double totalLoglik = topGtLoglik;
    ciCandidates.emplace(xFrom, xTo, totalLoglik);

    const int kMaxIntervalWidth = 750;
    double likelihoodRatio = 1;
    while (likelihoodRatio >= 0.01 && xTo - xFrom <= kMaxIntervalWidth)
    {
        double llShiftLeft0 = (genotyper->*likelihood)(xFrom - 1, yFrom - 1);
        double llShiftLeft1 = (genotyper->*likelihood)(xFrom - 1, yFrom);
        double llShiftLeft2 = (genotyper->*likelihood)(xFrom - 1, yFrom + 1);
        double llShiftLeft = std::max(std::max(llShiftLeft0, llShiftLeft1), llShiftLeft2);

        double llShiftRight0 = (genotyper->*likelihood)(xTo + 1, yTo + 1);
        double llShiftRight1 = (genotyper->*likelihood)(xTo + 1, yTo);
        double llShiftRight2 = (genotyper->*likelihood)(xTo + 1, yTo - 1);
        double llShiftRight = std::max(std::max(llShiftRight0, llShiftRight1), llShiftRight2);

        double gtLoglik;
        if (llShiftLeft >= llShiftRight)
        {
            xFrom -= 1;
            if (llShiftLeft0 > llShiftLeft1 && llShiftLeft0 > llShiftLeft2)
            {
                yFrom -= 1;
            }
            else if (llShiftLeft2 > llShiftLeft0 && llShiftLeft2 > llShiftLeft1)
            {
                yFrom += 1;
            }

            gtLoglik = (genotyper->*likelihood)(xFrom, yFrom);
        }
        else
        {
            xTo += 1;
            if (llShiftRight0 > llShiftRight1 && llShiftRight0 > llShiftRight2)
            {
                yTo += 1;
            }
            else if (llShiftRight2 > llShiftRight0 && llShiftRight2 > llShiftRight1)
            {
                yTo -= 1;
            }

            gtLoglik = (genotyper->*likelihood)(xTo, yTo);
        }

        totalLoglik = getLogSum(totalLoglik, gtLoglik);
        ciCandidates.emplace(xFrom, xTo, totalLoglik);
        likelihoodRatio = std::exp(gtLoglik - topGtLoglik);
    }

    auto ciCandidate = ciCandidates.top();
    while (ciCandidates.size() > 1)
    {
        ciCandidates.pop();
        auto nextCiCandidate = ciCandidates.top();
        const double ciProbability = std::exp(nextCiCandidate.loglik - totalLoglik);
        if (ciProbability >= 0.95)
        {
            ciCandidate = nextCiCandidate;
        }
        else
        {
            break;
        }
    }

    return { ciCandidate.startSize, ciCandidate.endSize };
}

RepeatGenotype TwoAlleleGenotyper::genotype(const unordered_set<int>& alleleSizeCandidates)
{
    RepeatGenotype gt = getMostLikelyGenotype(alleleSizeCandidates);
    const int bestShortSize = gt.shortAlleleSizeInUnits();
    const int bestLongSize = gt.longAlleleSizeInUnits();
    Ci shortStrCi = getCiAlongX(bestShortSize, bestLongSize, this, &TwoAlleleGenotyper::getShortAndLongAlleleLoglik);
    Ci longStrCi = getCiAlongX(bestLongSize, bestShortSize, this, &TwoAlleleGenotyper::getLongAndShortAlleleLoglik);

    gt.setShortAlleleSizeInUnitsCi(shortStrCi.begin, shortStrCi.end);
    gt.setLongAlleleSizeInUnitsCi(longStrCi.begin, longStrCi.end);

    return gt;
}

RepeatGenotype TwoAlleleGenotyper::getMostLikelyGenotype(const unordered_set<int>& alleleSizeCandidates)
{
    double maxGtLoglik = std::numeric_limits<double>::lowest();
    int bestShortAlleleSize = 0;
    int bestLongAlleleSize = 0;

    for (int shortAlleleSize : alleleSizeCandidates)
    {
        for (int longAlleleSize : alleleSizeCandidates)
        {
            if (shortAlleleSize > longAlleleSize)
            {
                continue;
            }

            const double gtLoglik = getShortAndLongAlleleLoglik(shortAlleleSize, longAlleleSize);

            if (maxGtLoglik < gtLoglik)
            {
                maxGtLoglik = gtLoglik;
                bestShortAlleleSize = shortAlleleSize;
                bestLongAlleleSize = longAlleleSize;
            }
        }
    }

    return RepeatGenotype(motifLen_, { bestShortAlleleSize, bestLongAlleleSize });
}

double TwoAlleleGenotyper::getShortAndLongAlleleLoglik(int shortAlleleSize, int longAlleleSize)
{
    if (shortAlleleSize < 0 || longAlleleSize < 0 || shortAlleleSize > longAlleleSize)
    {
        return std::numeric_limits<double>::lowest();
    }

    const int shortAlleleLen = shortAlleleSize * motifLen_ + fragLen_ + 1;
    const int longAlleleLen = longAlleleSize * motifLen_ + fragLen_ + 1;
    const double shortAlleleFrac = static_cast<double>(shortAlleleLen) / (shortAlleleLen + longAlleleLen);
    double genotypeLoglik = 0;

    for (int fragIndex = 0; fragIndex != fragLogliks_.numFrags(); ++fragIndex)
    {
        double fragLoglikForShortAllele = fragLogliks_.getLoglik(fragIndex, shortAlleleSize);
        double fragLoglikForLongAllele = fragLogliks_.getLoglik(fragIndex, longAlleleSize);

        const double shortAlleleTerm = std::log(shortAlleleFrac) + fragLoglikForShortAllele;
        const double longAlleleTerm = std::log(1.0 - shortAlleleFrac) + fragLoglikForLongAllele;

        const double loglikGivenRightmap = getLogSum(shortAlleleTerm, longAlleleTerm);
        const double loglikGivenMismap = topFragLogliks_[fragIndex];
        const double mismapPrior = std::log(0.001);
        const double rightmapPrior = std::log(1.0 - 0.001);
        genotypeLoglik += getLogSum(mismapPrior + loglikGivenMismap, rightmapPrior + loglikGivenRightmap);
    }

    return genotypeLoglik;
}

double TwoAlleleGenotyper::getLongAndShortAlleleLoglik(int longAlleleSize, int shortAlleleSize)
{
    return getShortAndLongAlleleLoglik(shortAlleleSize, longAlleleSize);
}

}
}
