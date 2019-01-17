//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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
