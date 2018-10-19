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

#include "reads/read.h"
#include "region_analysis/RegionAnalyzer.hh"

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

using IntervalWithLocusTypeAndAnalyzer = Interval<std::size_t, LocusTypeAndAnalyzer>;
using AnalyzerIntervalTree = IntervalTree<std::size_t, LocusTypeAndAnalyzer>;
using AnalyzerIntervalTrees = std::unordered_map<std::string, AnalyzerIntervalTree>;

class LocationBasedAnalyzerFinder
{
public:
    LocationBasedAnalyzerFinder(std::vector<std::unique_ptr<RegionAnalyzer>>& locusAnalyzers, int searchRadius);
    boost::optional<LocusTypeAndAnalyzer>
    query(const std::string& readChrom, int32_t readPosition, const std::string& mateChrom, int32_t matePosition);

private:
    boost::optional<LocusTypeAndAnalyzer>
    tryGettingLocusAnalyzer(const std::string& readChrom, int32_t readPosition) const;

    AnalyzerIntervalTrees intervalTrees_;
};
