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

#include <map>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "thirdparty/intervaltree/IntervalTree.h"

#include "reads/Read.hh"
#include "region_analysis/RegionAnalyzer.hh"

namespace ehunter
{

enum class LocusType
{
    kOfftargetLocus,
    kTargetLocus
};

struct LocusTypeAndAnalyzer
{
    LocusTypeAndAnalyzer(LocusType locusType, RegionAnalyzer* locusAnalyzerPtr)
        : locusType(locusType)
        , locusAnalyzerPtr(locusAnalyzerPtr)
    {
    }
    LocusType locusType;
    RegionAnalyzer* locusAnalyzerPtr;
};

using IntervalWithLocusTypeAndAnalyzer = ehunter::Interval<std::size_t, LocusTypeAndAnalyzer>;
using AnalyzerIntervalTree = ehunter::IntervalTree<std::size_t, LocusTypeAndAnalyzer>;
using AnalyzerIntervalTrees = std::unordered_map<int32_t, AnalyzerIntervalTree>;

class LocationBasedAnalyzerFinder
{
public:
    LocationBasedAnalyzerFinder(std::vector<std::unique_ptr<RegionAnalyzer>>& locusAnalyzers);
    boost::optional<LocusTypeAndAnalyzer>
    query(int32_t readContigId, int64_t readPosition, int32_t mateContigId, int64_t matePosition);

private:
    boost::optional<LocusTypeAndAnalyzer> tryGettingLocusAnalyzer(int32_t contigIndex, int64_t position) const;

    AnalyzerIntervalTrees intervalTrees_;
};

}
