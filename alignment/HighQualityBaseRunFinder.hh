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

#pragma once

#include <iostream>
#include <string>
#include <utility>

// Performs search for stretches of high-quality bases
class HighQualityBaseRunFinder
{
public:
    using StringIterPair = std::pair<std::string::const_iterator, std::string::const_iterator>;

    /**
     * @param windowSize: size of the window used for scanning an input sequence
     * @param minHighQualityBasesInGoodWindow: a window containing this many high-quality bases or more is deemed "good"
     * @param minLengthOfFullSizeRun: number bases in a run that is sufficent for the run to be reported
     */
    HighQualityBaseRunFinder(int windowSize, int minHighQualityBasesInGoodWindow, int minLengthOfFullSizeRun)
        : windowSize_(windowSize)
        , minHighQualityBasesInGoodWindow_(minHighQualityBasesInGoodWindow)
        , minLengthOfFullSizeRun_(minLengthOfFullSizeRun)
    {
    }

    /**
     * Searches for the first sufficiently-long run of high-quality bases
     *
     * @param query: any query sequence
     * @return pair of iterators delineating the run
     */
    StringIterPair find(const std::string& query) const;

private:
    std::string::const_iterator getStartOfNextBadWindow(
        std::string::const_iterator windowStartIter, std::string::const_iterator windowEndIter) const;
    bool isStartOfGoodWindow(std::string::const_iterator windowStartIter) const;

    const int windowSize_;
    const int minHighQualityBasesInGoodWindow_;
    const int minLengthOfFullSizeRun_;
};
