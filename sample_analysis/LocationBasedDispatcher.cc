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

#include "sample_analysis/LocationBasedDispatcher.hh"

using std::string;
using std::vector;

LocationBasedDispatcher::LocationBasedDispatcher(
    std::vector<std::unique_ptr<RegionAnalyzer>>& locusAnalyzers, int searchRadius)
    : locationBasedAnalyzerFinder_(locusAnalyzers, searchRadius)
{
}

void LocationBasedDispatcher::dispatch(
    const std::string& readChrom, int32_t readPosition, const std::string& mateChrom, int32_t matePosition,
    reads::Read read)
{
    // Check if the read is in the hash and act accordingly
    const auto mateIterator = unpairedReads_.find(read.fragmentId());
    if (mateIterator == unpairedReads_.end())
    {
        unpairedReads_.emplace(std::make_pair(read.fragmentId(), std::move(read)));
        return;
    }

    reads::Read& mate = mateIterator->second;
    auto optionalRegionTypeAndAnalyzer
        = locationBasedAnalyzerFinder_.query(readChrom, readPosition, mateChrom, matePosition);

    if (optionalRegionTypeAndAnalyzer && optionalRegionTypeAndAnalyzer->locusType == LocusType::kTargetLocus)
    {
        string fragmentId = mate.fragmentId();
        optionalRegionTypeAndAnalyzer->locusAnalyzerPtr->processMates(std::move(read), std::move(mate));
        unpairedReads_.erase(fragmentId);
    }
    else
    {
        unpairedReads_.erase(mate.fragmentId());
    }
}
