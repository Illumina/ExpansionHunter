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
#pragma once

#include "common/Common.hh"
#include <vector>

namespace ehunter
{

struct RegionDepthInfo
{
    RegionDepthInfo(double gc, double depth)
        : gc(gc)
        , depth(depth)
    {
    }
    double gc;
    double depth;
};

class DepthNormalizer
{
public:
    // given the GC and depth values of the normalization regions
    explicit DepthNormalizer(std::vector<RegionDepthInfo> normalizationRegions);

    // correct the original depth value by the difference between expected depth based on GC and the median of all
    // regions
    double correctDepth(double regionGC, double regionDepth, bool correctByGC) const;

    // exposed for unit testing
    const std::vector<double>& fittedDepths() const { return fittedDepths_; }

private:
    // scale factor for determining the amount of correction for the original depth
    double gcCorrectionScaleFactor_ = 0.9;
    // parameters for Lowess
    double smoothingSpan_ = 2.0 / 3.0;
    int numberIteration_ = 3;
    double deltaSkippingParameter_ = 0;

    // if too few regions for normalization, do not perform depth normalizatoin
    int minimumNumberOfNormalizatonRegions_ = 5;
    // if two flanking GC values are really close,
    // do not do interpolation and just take the left depth value instead
    double minimumGCRangeForInterpolation_ = 1e-4;

    double medianFittedDepth_ = 0;
    double medianDepth_ = 0;
    std::vector<double> fittedDepths_;
    std::vector<double> fittedGCs_;
    std::vector<double> robustnessWeights_;
    std::vector<double> residuals_;
};

double getMedian(std::vector<double> values);
}
