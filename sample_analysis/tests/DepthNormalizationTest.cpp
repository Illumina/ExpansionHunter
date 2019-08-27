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
#include "gtest/gtest.h"
#include <iostream>
#include <numeric>
#include <vector>

#include "common/Common.hh"

using namespace ehunter;
using std::vector;

TEST(DepthNormalizationByGC, TestMedian)
{
    vector<double> depthValuesOdd{ 18, 2, 15, 6, 10 };
    EXPECT_EQ(10, getMedian(depthValuesOdd));

    vector<double> depthValuesEven{ 2, 15, 6, 10 };
    EXPECT_EQ(8, getMedian(depthValuesEven));
}

TEST(DepthNormalizer, TestInitialization)
{
    vector<RegionDepthInfo> normalizationRegions
        = { RegionDepthInfo(0.51, 0.7), RegionDepthInfo(0.42, 0.88), RegionDepthInfo(0.49, 0.99),
            RegionDepthInfo(0.2, 1.05), RegionDepthInfo(0.4, 0.8) };
    DepthNormalizer normalizer(normalizationRegions);

    vector<double> expectedFittedDepths{ 1.193, 0.909, 1, 1.125, 0.795 };

    for (int index = 0; index != static_cast<int>(expectedFittedDepths.size()); index++)
    {
        EXPECT_NEAR(expectedFittedDepths[index], normalizer.fittedDepths()[index], 1e-3);
    }
}

TEST(DepthNormalizationByGC, TestCorrectByGC)
{

    {
        DepthNormalizer normalizer({ RegionDepthInfo(0.51, 0.7), RegionDepthInfo(0.42, 0.88) });
        double depthValue = 0.8;
        double gcValueExact = 0.4;
        EXPECT_NEAR(0.8, normalizer.correctDepth(gcValueExact, depthValue, true), 1e-3);
    }

    {
        DepthNormalizer normalizer({ RegionDepthInfo(0.51, 0.7), RegionDepthInfo(0.42, 0.88),
                                     RegionDepthInfo(0.49, 0.99), RegionDepthInfo(0.2, 1.05),
                                     RegionDepthInfo(0.4, 0.8) });

        double depthValue = 0.8;

        double gcValueExact = 0.4;
        EXPECT_NEAR(0.909, normalizer.correctDepth(gcValueExact, depthValue, false), 1e-3);
        EXPECT_NEAR(0.983, normalizer.correctDepth(gcValueExact, depthValue, true), 1e-3);

        double gcValueLow = 0.1;
        EXPECT_NEAR(0.751, normalizer.correctDepth(gcValueLow, depthValue, true), 1e-3);

        double gcValueHigh = 0.6;
        EXPECT_NEAR(1.076, normalizer.correctDepth(gcValueHigh, depthValue, true), 1e-3);

        double gcValueInexact = 0.455;
        EXPECT_NEAR(0.858, normalizer.correctDepth(gcValueInexact, depthValue, true), 1e-3);
    }
}
