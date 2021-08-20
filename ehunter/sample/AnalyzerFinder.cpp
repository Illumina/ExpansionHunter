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

#include "sample/AnalyzerFinder.hh"

using boost::optional;
using ehunter::locus::LocusAnalyzer;
using ehunter::locus::RegionType;
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

/// Remove items from \p bundles which refer to a LocusAnalyzer already found in \p commonBundles
///
void filterOutCommonBundles(const vector<AnalyzerBundle>& commonBundles, vector<AnalyzerBundle>& bundles)
{
    if (commonBundles.size() == bundles.size())
    {
        bundles.clear();
    }
    else
    {
        vector<AnalyzerBundle> remainingBundles;
        for (auto& bundle : bundles)
        {
            bool bundleAssigned(false);
            for (const auto& commonBundle : commonBundles)
            {
                if (bundle.locusIndex == commonBundle.locusIndex)
                {
                    bundleAssigned = true;
                    break;
                }
            }
            if (not bundleAssigned)
            {
                remainingBundles.emplace_back(bundle);
            }
        }
        bundles = remainingBundles;
    }
}

/// \param[in,out] readBundles Used to find common bundles, updated to contain any remaining bundles on return
///
/// \param[in,out] mateBundles Used to find common bundles, updated to contain any remaining bundles on return
///
vector<AnalyzerBundle> coalesceCommonBundles(vector<AnalyzerBundle>& readBundles, vector<AnalyzerBundle>& mateBundles)
{
    vector<AnalyzerBundle> commonBundles;

    for (const auto& readBundle : readBundles)
    {
        for (const auto& mateBundle : mateBundles)
        {
            if (readBundle.locusIndex == mateBundle.locusIndex)
            {
                commonBundles.emplace_back(readBundle);
                commonBundles.back().regionType = coalesceRegionTypes(readBundle.regionType, mateBundle.regionType);
                break;
            }
        }
    }

    if (not commonBundles.empty())
    {
        filterOutCommonBundles(commonBundles, readBundles);
        filterOutCommonBundles(commonBundles, mateBundles);
    }

    return commonBundles;
}

/// We ignore nearby pairs where one mate is inside and one mate is outside of the offtarget region
///
/// \param[in,out] bundles Coalesced bundles are appended to this structure
///
void coalesceBundlesForNearbyMates(
    const vector<AnalyzerBundle>& readBundles, const vector<AnalyzerBundle>& mateBundles,
    vector<AnalyzerBundle>& bundles)
{
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
}

/// \param[in,out] bundles Coalesced bundles are appended to this structure
///
void coalesceBundlesForFarawayMates(
    const vector<AnalyzerBundle>& readBundles, const vector<AnalyzerBundle>& mateBundles,
    vector<AnalyzerBundle>& bundles)
{
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
}
}

void processAnalyzerBundleReadPair(
    locus::LocusAnalyzer& locusAnalyzer, locus::RegionType regionType, AnalyzerInputType inputType, Read& read,
    Read& mate, graphtools::AlignerSelector& alignerSelector)
{

    switch (inputType)
    {
    case AnalyzerInputType::kBothReads:
        locusAnalyzer.processMates(read, &mate, regionType, alignerSelector);
        break;
    case AnalyzerInputType::kReadOnly:
        locusAnalyzer.processMates(read, nullptr, regionType, alignerSelector);
        break;
    case AnalyzerInputType::kMateOnly:
        locusAnalyzer.processMates(mate, nullptr, regionType, alignerSelector);
        break;
    }
}

AnalyzerFinder::AnalyzerFinder(vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers)
{
    using IntervalWithLocusTypeAndAnalyzer = Interval<std::size_t, AnalyzerBundle>;

    unordered_map<int32_t, vector<IntervalWithLocusTypeAndAnalyzer>> contigToIntervals;

    const unsigned locusAnalzerCount(locusAnalyzers.size());
    for (unsigned locusAnalyzerIndex(0); locusAnalyzerIndex < locusAnalzerCount; ++locusAnalyzerIndex)
    {
        const auto& locusAnalyzer(locusAnalyzers[locusAnalyzerIndex]);
        const LocusSpecification& locusSpec = locusAnalyzer->locusSpec();
        for (const auto& region : locusSpec.targetReadExtractionRegions())
        {
            AnalyzerBundle bundle(RegionType::kTarget, locusAnalyzerIndex);
            contigToIntervals[region.contigIndex()].emplace_back(region.start(), region.end(), bundle);
        }

        for (const auto& region : locusSpec.offtargetReadExtractionRegions())
        {
            AnalyzerBundle bundle(RegionType::kOfftarget, locusAnalyzerIndex);
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
    vector<AnalyzerBundle> bundles = coalesceCommonBundles(readAnalyzerBundles, mateAnalyzerBundles);

    if ((not readAnalyzerBundles.empty()) or (not mateAnalyzerBundles.empty()))
    {
        if (areMatesNearby(readContigId, readStart, mateContigId, mateStart))
        {
            coalesceBundlesForNearbyMates(readAnalyzerBundles, mateAnalyzerBundles, bundles);
        }
        else
        {
            coalesceBundlesForFarawayMates(readAnalyzerBundles, mateAnalyzerBundles, bundles);
        }
    }
    return bundles;
}

}
