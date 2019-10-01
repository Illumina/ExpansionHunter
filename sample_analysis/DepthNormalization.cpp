//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Xiao Chen <xchen2@illumina.com>
//         Egor Dolzhenko <edolzhenko@illumina.com>
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
#include "sample_analysis/DepthNormalization.hh"
#include "stats/LowessRegression.hh"
#include <algorithm>
#include <cassert>
#include <iostream>
using std::vector;

namespace ehunter
{

double getMedian(std::vector<double> values)
{
    int vectorSize = static_cast<int>(values.size());
    std::sort(values.begin(), values.end());
    int medianPosition1 = vectorSize / 2;
    int medianPosition2 = vectorSize - vectorSize / 2 - 1;
    double median = 0.5 * (values[medianPosition1] + values[medianPosition2]);

    return median;
}

static vector<double> getDepths(const vector<RegionDepthInfo>& regionInfos)
{
    vector<double> depths;
    depths.reserve(regionInfos.size());
    for (const auto& regionInfo : regionInfos)
    {
        depths.push_back(regionInfo.depth);
    }
    return depths;
}

static vector<double> getGCs(const vector<RegionDepthInfo>& regionInfos)
{
    vector<double> gsc;
    gsc.reserve(regionInfos.size());
    for (const auto& regionInfo : regionInfos)
    {
        gsc.push_back(regionInfo.gc);
    }
    return gsc;
}

DepthNormalizer::DepthNormalizer(vector<RegionDepthInfo> normalizationRegions)
{
    // Normalize by median depth
    vector<double> depths = getDepths(normalizationRegions);
    medianDepth_ = getMedian(depths);
    for (auto& regionInfo : normalizationRegions)
    {
        regionInfo.depth = regionInfo.depth / medianDepth_;
    }

    // Sort based on GC
    std::sort(
        normalizationRegions.begin(), normalizationRegions.end(),
        [](const RegionDepthInfo& r1, const RegionDepthInfo& r2) { return r1.gc < r2.gc; });

    fittedGCs_ = getGCs(normalizationRegions);
    fittedDepths_.resize(normalizationRegions.size());
    robustnessWeights_.resize(normalizationRegions.size());
    residuals_.resize(normalizationRegions.size());

    LowessRegression lowessRegresser(smoothingSpan_, deltaSkippingParameter_, numberIteration_);
    lowessRegresser.regression(
        fittedGCs_, getDepths(normalizationRegions), fittedDepths_, robustnessWeights_, residuals_);
    medianFittedDepth_ = getMedian(fittedDepths_);
}

double DepthNormalizer::correctDepth(double regionGC, double regionDepth, bool correctByGC) const
{
    // no correction
    if (static_cast<int>(fittedDepths_.size()) < minimumNumberOfNormalizatonRegions_)
    {
        return regionDepth;
    }
    // normalize first by sample median
    regionDepth /= medianDepth_;
    if (!correctByGC)
    {
        return regionDepth;
    }
    double expectedDepthForGivenGC;
    auto findItem = std::find(fittedGCs_.begin(), fittedGCs_.end(), regionGC);
    // found the exact GC value in normalization regions
    if (findItem != fittedGCs_.end())
    {
        auto itemIndex = distance(fittedGCs_.begin(), findItem);
        expectedDepthForGivenGC = fittedDepths_[itemIndex];
    }
    else
    {
        auto upperBound = std::upper_bound(fittedGCs_.begin(), fittedGCs_.end(), regionGC);
        // GC lower than all normalization regions
        if (upperBound == fittedGCs_.begin())
        {
            expectedDepthForGivenGC = fittedDepths_[0];
        }
        // GC higher than all normalization regions
        else if (upperBound == fittedGCs_.end())
        {
            expectedDepthForGivenGC = fittedDepths_[fittedDepths_.size() - 1];
        }
        // interpolate using the two flanking GC values
        else
        {
            auto upperBoundIndex = distance(fittedGCs_.begin(), upperBound);
            auto lowerBoundIndex = upperBoundIndex - 1;
            double gcRange = fittedGCs_[upperBoundIndex] - fittedGCs_[lowerBoundIndex];
            if (gcRange < minimumGCRangeForInterpolation_)
            {
                expectedDepthForGivenGC = fittedGCs_[lowerBoundIndex];
            }
            else
            {
                double ratio = (regionGC - fittedGCs_[lowerBoundIndex]) / gcRange;
                expectedDepthForGivenGC
                    = ratio * fittedDepths_[upperBoundIndex] + (1 - ratio) * fittedDepths_[lowerBoundIndex];
            }
        }
    }

    // correct depth by the difference between GC-expected depth and sample median, modified by a scale factor
    double depthScaleFactor = gcCorrectionScaleFactor_ * std::min(regionDepth, 2.0);
    double gcCorrectedDepth = regionDepth + depthScaleFactor * (medianFittedDepth_ - expectedDepthForGivenGC);
    return gcCorrectedDepth;
}
}