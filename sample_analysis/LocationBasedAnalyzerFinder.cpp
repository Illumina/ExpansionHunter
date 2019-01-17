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

#include "sample_analysis/LocationBasedAnalyzerFinder.hh"

#include <stdexcept>
#include <unordered_map>

using boost::optional;
using std::map;
using std::string;
using std::unordered_map;
using std::vector;

namespace ehunter
{

LocationBasedAnalyzerFinder::LocationBasedAnalyzerFinder(
    vector<std::unique_ptr<RegionAnalyzer>>& locusAnalyzers, int searchRadius)
{
    unordered_map<string, vector<IntervalWithLocusTypeAndAnalyzer>> contigToIntervals;
    for (auto& locusAnalyzer : locusAnalyzers)
    {
        const LocusSpecification& locusSpec = locusAnalyzer->regionSpec();
        for (auto& refRegion : locusSpec.referenceLoci())
        {
            auto targetRegion = refRegion.extend(searchRadius);
            LocusTypeAndAnalyzer payload(LocusType::kTargetLocus, locusAnalyzer.get());
            contigToIntervals[targetRegion.chrom()].emplace_back(targetRegion.start(), targetRegion.end(), payload);
        }
        for (const auto& offtargetLocus : locusSpec.offtargetLoci())
        {
            LocusTypeAndAnalyzer payload(LocusType::kOfftargetLocus, locusAnalyzer.get());
            contigToIntervals[offtargetLocus.chrom()].emplace_back(
                offtargetLocus.start(), offtargetLocus.end(), payload);
        }
    }

    for (auto& contigAndIntervals : contigToIntervals)
    {
        const string& contig = contigAndIntervals.first;
        auto intervals = contigAndIntervals.second;
        intervalTrees_.emplace(std::make_pair(contig, AnalyzerIntervalTree(std::move(intervals))));
    }
}

optional<LocusTypeAndAnalyzer> LocationBasedAnalyzerFinder::query(
    const string& readChrom, int32_t readPosition, const string& mateChrom, int32_t matePosition)
{
    auto optionalReadLocusTypeAndAnalyzer = tryGettingLocusAnalyzer(readChrom, readPosition);
    auto optionalMateLocusTypeAndAnalyzer = tryGettingLocusAnalyzer(mateChrom, matePosition);

    if (optionalReadLocusTypeAndAnalyzer && optionalReadLocusTypeAndAnalyzer->locusType == LocusType::kTargetLocus)
    {
        return optionalReadLocusTypeAndAnalyzer;
    }

    if (optionalMateLocusTypeAndAnalyzer && optionalMateLocusTypeAndAnalyzer->locusType == LocusType::kTargetLocus)
    {
        return optionalMateLocusTypeAndAnalyzer;
    }

    if (optionalReadLocusTypeAndAnalyzer)
    {
        return optionalReadLocusTypeAndAnalyzer;
    }

    if (optionalMateLocusTypeAndAnalyzer)
    {
        return optionalMateLocusTypeAndAnalyzer;
    }

    return optional<LocusTypeAndAnalyzer>();
}

optional<LocusTypeAndAnalyzer>
LocationBasedAnalyzerFinder::tryGettingLocusAnalyzer(const string& readChrom, int32_t readPosition) const
{
    const auto contigTreeIterator = intervalTrees_.find(readChrom);

    if (contigTreeIterator == intervalTrees_.end())
    {
        return optional<LocusTypeAndAnalyzer>();
    }

    const vector<IntervalWithLocusTypeAndAnalyzer>& relevantRegionAnalyzers
        = contigTreeIterator->second.findOverlapping(readPosition, readPosition + 1);

    if (relevantRegionAnalyzers.size() > 1)
    {
        throw std::logic_error("Repeat catalog must contain non-overlapping regions");
    }

    if (!relevantRegionAnalyzers.empty())
    {
        return relevantRegionAnalyzers.front().value;
    }
    else
    {
        return optional<LocusTypeAndAnalyzer>();
    }
}

}
