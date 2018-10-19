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

using std::string;

HighQualityBaseRunFinder::StringIterPair HighQualityBaseRunFinder::find(const string& query) const
{
    auto endOfWindowStarts = query.end() - windowSize_ + 1;
    auto currentWindowStart = query.begin();
    auto startOfNextBadWindow = currentWindowStart;

    while (currentWindowStart != endOfWindowStarts)
    {
        startOfNextBadWindow = getStartOfNextBadWindow(currentWindowStart, endOfWindowStarts);

        const int runLength = startOfNextBadWindow - currentWindowStart;
        if (runLength >= minLengthOfFullSizeRun_ || startOfNextBadWindow == endOfWindowStarts)
        {
            break;
        }
        else if (startOfNextBadWindow == currentWindowStart)
        {
            ++currentWindowStart;
        }
        else
        {
            currentWindowStart = startOfNextBadWindow;
        }
    }

    if (startOfNextBadWindow == endOfWindowStarts)
    {
        startOfNextBadWindow = query.end();
    }

    const int runLength = startOfNextBadWindow - currentWindowStart;
    if (runLength < minLengthOfFullSizeRun_)
    {
        currentWindowStart = query.end();
        startOfNextBadWindow = query.end();
    }

    return std::make_pair(currentWindowStart, startOfNextBadWindow);
}

string::const_iterator HighQualityBaseRunFinder::getStartOfNextBadWindow(
    string::const_iterator currentWindowStart, string::const_iterator windowEndIter) const
{
    while (currentWindowStart != windowEndIter)
    {
        if (!isStartOfGoodWindow(currentWindowStart))
        {
            return currentWindowStart;
        }
        ++currentWindowStart;
    }

    return windowEndIter;
}

bool HighQualityBaseRunFinder::isStartOfGoodWindow(string::const_iterator currentWindowStart) const
{
    int numHighQualityBases = 0;

    auto windowEndIter = currentWindowStart + windowSize_;
    for (auto windowIter = currentWindowStart; windowIter != windowEndIter; ++windowIter)
    {
        if (isupper(*windowIter))
        {
            if (++numHighQualityBases == minHighQualityBasesInGoodWindow_)
            {
                return true;
            }
        }
    }

    return false;
}
