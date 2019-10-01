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

#include "stats/LowessRegression.hh"
#include "gtest/gtest.h"
#include <numeric>
#include <vector>

#include "common/Common.hh"

using namespace ehunter;

TEST(Lowess, TestLowess)
{
    std::vector<double> testYValues{ 18, 2, 15, 6, 10, 4, 16, 11, 7, 3, 14, 17, 20, 12, 9, 13, 1, 8, 5, 19 };
    std::vector<double> testXValues{ 1, 2, 3, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 8, 10, 12, 14, 50 };
    std::vector<double> expectedResult{
        13.659, 11.145, 8.701,  9.722,  10.000, 11.300, 11.300, 11.300, 11.300, 11.300,
        11.300, 11.300, 11.300, 11.300, 11.300, 13.000, 6.440,  5.596,  5.456,  18.998
    };

    LowessRegression lowessRegresser(0.25, 0, 0);
    int vectorSize = testXValues.size();
    std::vector<double> fittedYValues(vectorSize);
    std::vector<double> robustnessWeights(vectorSize);
    std::vector<double> fitResiduals(vectorSize);

    lowessRegresser.regression(testXValues, testYValues, fittedYValues, robustnessWeights, fitResiduals);

    for (int i = 0; i != (int)expectedResult.size(); i++)
    {
        EXPECT_NEAR(expectedResult[i], fittedYValues[i], 1e-3);
    }
}