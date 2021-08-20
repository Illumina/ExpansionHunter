//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "alignment/HighQualityBaseRunFinder.hh"

#include <cctype>
#include <limits>

using std::string;

namespace ehunter
{

template <typename Iter> static int countGoodBases(Iter begin, Iter end)
{
    int numGoodBases = 0;
    for (; begin != end; ++begin)
    {
        if (isupper(*begin))
        {
            ++numGoodBases;
        }
    }

    return numGoodBases;
}

template <typename Iter> static double calculateBaseRunProb(double goodBaseProb, Iter begin, Iter end)
{
    const int runLength = end - begin;
    const int numGoodBasesInRun = countGoodBases(begin, end);
    const int numBadBasesInRun = runLength - numGoodBasesInRun;
    return numGoodBasesInRun * goodBaseProb + numBadBasesInRun * (1.0 - goodBaseProb);
}

template <typename Iter>
static double calculateRunProbability(
    double probOfGoodBaseInFirstRun, double probOfGoodBaseInSecondRun, Iter start, Iter changePoint, Iter end)
{
    const double firstRunProb = calculateBaseRunProb(probOfGoodBaseInFirstRun, start, changePoint);
    const double secondRunProb = calculateBaseRunProb(probOfGoodBaseInSecondRun, changePoint, end);

    return firstRunProb + secondRunProb;
}

template <typename Iter>
static Iter findTopChangePoint(double probOfGoodBaseInFirstRun, double probOfGoodBaseInSecondRun, Iter start, Iter end)
{
    double topRunProb = 0;
    auto topChangePoint = start;

    for (auto changePoint = start; changePoint != end; ++changePoint)
    {
        const double currentRunProb
            = calculateRunProbability(probOfGoodBaseInFirstRun, probOfGoodBaseInSecondRun, start, changePoint, end);

        if (topRunProb < currentRunProb)
        {
            topRunProb = currentRunProb;
            topChangePoint = changePoint;
        }
    }

    return topChangePoint;
}

StringIterPair
findHighQualityBaseRun(const string& query, double probOfGoodBaseInBadRun, double probOfGoodBaseInGoodRun)
{
    const auto middleBaseIter = query.begin() + query.length() / 2;
    const auto startOfGoodRun
        = findTopChangePoint(probOfGoodBaseInBadRun, probOfGoodBaseInGoodRun, query.begin(), middleBaseIter);

    const auto middleBaseReverseIter = query.rbegin() + query.length() / 2;
    const auto topEndingChangePointReverseIter
        = findTopChangePoint(probOfGoodBaseInBadRun, probOfGoodBaseInGoodRun, query.rbegin(), middleBaseReverseIter);
    const auto endOfGoodRun = topEndingChangePointReverseIter.base();

    return std::make_pair(startOfGoodRun, endOfGoodRun);
}

}
