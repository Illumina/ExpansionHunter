//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
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

#include "sample_analysis/AnalyzerFinder.hh"

using boost::optional;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

namespace ehunter
{

namespace
{
    bool areMatesNearby(int32_t readContigId, int64_t readPosition, int32_t mateContigId, int64_t matePosition)
    {
        const int kMaxMateDistance = 1000;
        return ((readContigId == mateContigId) && (std::abs(readPosition - matePosition) < kMaxMateDistance));
    }

    inline RegionType coalesceRegionTypes(RegionType readRegionType, RegionType mateRegionType)
    {
        if (readRegionType == RegionType::kTarget || mateRegionType == RegionType::kTarget)
        {
            return RegionType::kTarget;
        }

        return RegionType::kOfftarget;
    }

    vector<AnalyzerBundle>
    coalesceCommonBundles(const vector<AnalyzerBundle>& readBundles, const vector<AnalyzerBundle>& mateBundles)
    {
        vector<AnalyzerBundle> commonBundles;

        for (const auto& readBundle : readBundles)
        {
            for (const auto& mateBundle : mateBundles)
            {
                if (readBundle.locusAnalyzerPtr == mateBundle.locusAnalyzerPtr)
                {
                    commonBundles.push_back(readBundle);
                    commonBundles.back().regionType = coalesceRegionTypes(readBundle.regionType, mateBundle.regionType);
                }
            }
        }

        return commonBundles;
    }

    // We ignore nearby pairs where one mate is inside and one mate is outside of the offtarget region
    vector<AnalyzerBundle>
    coalesceBundlesForNearbyMates(const vector<AnalyzerBundle>& readBundles, const vector<AnalyzerBundle>& mateBundles)
    {
        vector<AnalyzerBundle> bundles;

        for (const auto& bundle : readBundles)
        {
            if (bundle.regionType == RegionType::kTarget)
            {
                bundles.push_back(bundle);
                bundles.back().inputType = AnalyzerInputType::kReadOnly;
            }
        }

        for (const auto& bundle : mateBundles)
        {
            if (bundle.regionType == RegionType::kTarget)
            {
                bundles.push_back(bundle);
                bundles.back().inputType = AnalyzerInputType::kMateOnly;
            }
        }

        return bundles;
    }

    vector<AnalyzerBundle>
    coalesceBundlesForFarawayMates(const vector<AnalyzerBundle>& readBundles, const vector<AnalyzerBundle>& mateBundles)
    {
        vector<AnalyzerBundle> bundles;

        for (const auto& bundle : readBundles)
        {
            bundles.push_back(bundle);
            bundles.back().inputType = AnalyzerInputType::kBothReads;
        }

        for (const auto& bundle : mateBundles)
        {
            bundles.push_back(bundle);
            bundles.back().inputType = AnalyzerInputType::kBothReads;
        }

        return bundles;
    }
}

AnalyzerFinder::AnalyzerFinder(vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers)
{
    using IntervalWithLocusTypeAndAnalyzer = ehunter::Interval<std::size_t, AnalyzerBundle>;

    unordered_map<int32_t, vector<IntervalWithLocusTypeAndAnalyzer>> contigToIntervals;
    for (auto& locusAnalyzer : locusAnalyzers)
    {
        const LocusSpecification& locusSpec = locusAnalyzer->locusSpec();
        for (const auto& region : locusSpec.targetReadExtractionRegions())
        {
            AnalyzerBundle bundle(RegionType::kTarget, locusAnalyzer.get());
            contigToIntervals[region.contigIndex()].emplace_back(region.start(), region.end(), bundle);
        }

        for (const auto& region : locusSpec.offtargetReadExtractionRegions())
        {
            AnalyzerBundle bundle(RegionType::kOfftarget, locusAnalyzer.get());
            contigToIntervals[region.contigIndex()].emplace_back(region.start(), region.end(), bundle);
        }
    }

    for (auto& contigAndIntervals : contigToIntervals)
    {
        int32_t contigIndex = contigAndIntervals.first;
        auto intervals = contigAndIntervals.second;

        intervalTrees_.emplace(std::make_pair(contigIndex, AnalyzerIntervalTree(std::move(intervals))));
    }
}

vector<AnalyzerBundle> AnalyzerFinder::query(int32_t contigIndex, int64_t start, int64_t end) const
{
    const auto contigTreeIterator = intervalTrees_.find(contigIndex);
    if (contigTreeIterator == intervalTrees_.end())
    {
        return vector<AnalyzerBundle>();
    }

    const auto& intervalsWithBundles = contigTreeIterator->second.findOverlapping(start, end);

    vector<AnalyzerBundle> analyzerBundles;
    for (const auto& intervalWithBundle : intervalsWithBundles)
    {
        const bool isReadContainedInInterval = static_cast<int64_t>(intervalWithBundle.start) <= start
            && end <= static_cast<int64_t>(intervalWithBundle.stop);
        if (isReadContainedInInterval)
        {
            analyzerBundles.push_back(intervalWithBundle.value);
        }
    }

    return analyzerBundles;
}

vector<AnalyzerBundle> AnalyzerFinder::query(
    int32_t readContigId, int64_t readStart, int64_t readEnd, int32_t mateContigId, int64_t mateStart,
    int64_t mateEnd) const
{
    vector<AnalyzerBundle> readAnalyzerBundles = query(readContigId, readStart, readEnd);
    vector<AnalyzerBundle> mateAnalyzerBundles = query(mateContigId, mateStart, mateEnd);
    vector<AnalyzerBundle> commonBundles = coalesceCommonBundles(readAnalyzerBundles, mateAnalyzerBundles);

    if (!commonBundles.empty())
    {
        return commonBundles;
    }
    else if (areMatesNearby(readContigId, readStart, mateContigId, mateStart))
    {
        return coalesceBundlesForNearbyMates(readAnalyzerBundles, mateAnalyzerBundles);
    }
    else
    {
        return coalesceBundlesForFarawayMates(readAnalyzerBundles, mateAnalyzerBundles);
    }
}

}
