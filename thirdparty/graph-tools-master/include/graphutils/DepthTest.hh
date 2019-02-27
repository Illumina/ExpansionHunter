//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Sai Chen <schen6@illumina.com>,
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

#pragma once

#include <boost/math/distributions/normal.hpp>
#include <cstdint>

using boost::math::normal_distribution;

class DepthTest
{
public:
    /**
     * Set up depth test
     * @param expectedNumReads The mean number of reads when the coverage is as expected
     * @param stdDeviation The standard deviation of the number of reads
     * @param lowerSignificanceThreshold P value cutoff at the lower tail of the distribution
     * @param upperSignificanceThreshold P value cutoff at the upper tail of the distribution
     */
    DepthTest(
        int32_t expectedNumReads, double stdDeviation, double lowerSignificanceThreshold,
        double upperSignificanceThreshold);

    /**
     * Given observed number of reads, return true if pass depth test
     */
    bool testReadCount(int32_t observedNumReads);

private:
    int32_t expectedNumReads_;
    double stdDeviation_;
    double lowerSignificanceThreshold_;
    double upperSignificanceThreshold_;
    const normal_distribution<> coverageDistribution_;
};
