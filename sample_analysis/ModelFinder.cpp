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

#include "sample_analysis/ModelFinder.hh"

#include <memory>

#include "region_spec/LocusSpecification.hh"
#include "workflow/RegionModel.hh"

using boost::optional;
using std::shared_ptr;
using std::unique_ptr;
using std::unordered_map;
using std::unordered_set;
using std::vector;

namespace ehunter
{

namespace
{

    /*
    vector<AnalyzerBundle>
    coalesceCommonBundles(const vector<AnalyzerBundle>& readBundles, const vector<AnalyzerBundle>& mateBundles)
    {
        vector<AnalyzerBundle> commonBundles;

        for (const auto& readBundle : readBundles)
        {
            for (const auto& mateBundle : mateBundles)
            {
                if (readBundle.regionPtr == mateBundle.regionPtr)
                {
                    if (readBundle.regionPtr->type() == RegionModel::Type::kTarget)
                    {
                        commonBundles.push_back(readBundle);
                    }
                    else
                    {
                        commonBundles.push_back(mateBundle);
                    }
                }
            }
        }

        return commonBundles;
    } */

    // We ignore nearby pairs where one mate is inside and one mate is outside of the offtarget workflow
    /*
    vector<AnalyzerBundle>
    coalesceBundlesForNearbyMates(const vector<AnalyzerBundle>& readBundles, const vector<AnalyzerBundle>& mateBundles)
    {
        vector<AnalyzerBundle> bundles;

        for (const auto& bundle : readBundles)
        {
            if (bundle.regionPtr->type() == RegionModel::Type::kTarget)
            {
                bundles.push_back(bundle);
                bundles.back().inputType = AnalyzerInputType::kReadOnly;
            }
        }

        for (const auto& bundle : mateBundles)
        {
            if (bundle.regionPtr->type() == RegionModel::Type::kTarget)
            {
                bundles.push_back(bundle);
                bundles.back().inputType = AnalyzerInputType::kMateOnly;
            }
        }

        return bundles;
    } */

    /*
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
    } */
}

ModelFinder::ModelFinder(const vector<shared_ptr<RegionModel>>& models)
{
    using IntervalWithModel = ehunter::Interval<std::size_t, RegionModel*>;

    unordered_map<int, vector<IntervalWithModel>> contigToIntervals;
    for (auto& regionModelPtr : models)
    {
        for (const auto& genomicRegion : regionModelPtr->readExtractionRegions())
        {
            contigToIntervals[genomicRegion.contigIndex()].emplace_back(
                genomicRegion.start(), genomicRegion.end(), regionModelPtr.get());
        }
    }

    for (auto& contigAndIntervals : contigToIntervals)
    {
        int contigIndex = contigAndIntervals.first;
        auto intervals = contigAndIntervals.second;

        contigToModelTrees_.emplace(std::make_pair(contigIndex, ModelTree(std::move(intervals))));
    }
}

unordered_set<RegionModel*> ModelFinder::query(int contigId, int64_t start, int64_t end) const
{
    const auto iter = contigToModelTrees_.find(contigId);
    if (iter == contigToModelTrees_.end())
    {
        return unordered_set<RegionModel*>();
    }

    const auto& intervalsWithModels = iter->second.findOverlapping(start, end);

    unordered_set<RegionModel*> models;
    for (const auto& interval : intervalsWithModels)
    {
        // const bool fullContainment
        //    = static_cast<int64_t>(interval.start) <= start && end <= static_cast<int64_t>(interval.stop);
        // if (fullContainment)
        //{
        models.insert(interval.value);
        //}
    }

    return models;
}

/*
vector<AnalyzerBundle> ModelFinder::query(
    int readContigId, int64_t readStart, int64_t readEnd, int mateContigId, int64_t mateStart, int64_t mateEnd) const
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
} */

}
