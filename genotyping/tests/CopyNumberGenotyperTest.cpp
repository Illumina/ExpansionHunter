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

#include "genotyping/CopyNumberGenotyper.hh"
#include "gtest/gtest.h"
#include <numeric>
#include <vector>

#include "common/Common.hh"

using boost::optional;
using namespace ehunter;

TEST(CopyNumberGenotyping, ThrowsWithIllegalParameter)
{
    const int maxCopyNumber = 2;
    const double depthScaleFactor = 1;
    const double standardDeviationCN2 = 0.1;

    std::vector<double> priorFrequencies{ 0.1, 0.2, 0.7 };
    std::vector<double> tooFewMeanDepthValues{ 1, 2 };
    ASSERT_ANY_THROW(CopyNumberGenotyper(
        maxCopyNumber, depthScaleFactor, standardDeviationCN2, tooFewMeanDepthValues, priorFrequencies));

    std::vector<double> meanDepthValues{ 0, 1, 2 };
    std::vector<double> tooManyPriorFrequencies{ 0.1, 0.2, 0.3, 0.4 };
    ASSERT_ANY_THROW(CopyNumberGenotyper(
        maxCopyNumber, depthScaleFactor, standardDeviationCN2, meanDepthValues, tooManyPriorFrequencies));
}

TEST(CopyNumberGenotyping, RegularGenotype)
{
    int maxCopyNumber = 2;
    double depthScaleFactor = 1;
    double standardDeviationCN2 = 0.1;
    std::vector<double> priorFrequencies{ 0.1, 0.2, 0.7 };
    std::vector<double> meanDepthValues{ 0, 1, 2 };

    {
        CopyNumberGenotyper genotyper(
            maxCopyNumber, depthScaleFactor, standardDeviationCN2, meanDepthValues, priorFrequencies);
        EXPECT_EQ(2, *genotyper.genotype(2.05));
        EXPECT_EQ(1, *genotyper.genotype(0.95));
        EXPECT_EQ(0, *genotyper.genotype(0.01));
        EXPECT_EQ(boost::none, genotyper.genotype(3.05));
    }

    {
        // test a different scale factor
        depthScaleFactor = 2;
        CopyNumberGenotyper genotyper(
            maxCopyNumber, depthScaleFactor, standardDeviationCN2, meanDepthValues, priorFrequencies);
        EXPECT_EQ(1, *genotyper.genotype(2.05));
        EXPECT_EQ(boost::none, genotyper.genotype(0.95));
        EXPECT_EQ(0, *genotyper.genotype(0.01));
    }

    {
        // test different mean values
        depthScaleFactor = 1;
        meanDepthValues = { 0.5, 0.6, 2 };
        CopyNumberGenotyper genotyper(
            maxCopyNumber, depthScaleFactor, standardDeviationCN2, meanDepthValues, priorFrequencies);
        EXPECT_EQ(2, *genotyper.genotype(2.05));
        EXPECT_EQ(boost::none, genotyper.genotype(0.5));
        EXPECT_EQ(boost::none, genotyper.genotype(0.01));
    }

    {
        // more copy number states allowed
        maxCopyNumber = 3;
        priorFrequencies = { 0.1, 0.2, 0.3, 0.4 };
        meanDepthValues = { 0, 1, 2, 3 };
        CopyNumberGenotyper genotyper(
            maxCopyNumber, depthScaleFactor, standardDeviationCN2, meanDepthValues, priorFrequencies);
        EXPECT_EQ(2, *genotyper.genotype(2.05));
        EXPECT_EQ(3, *genotyper.genotype(3.05));
    }
}

TEST(CopyNumberGenotyping, BestGenotypeAndPosterior)
{
    const int maxCopyNumber = 2;
    const double depthScaleFactor = 1;
    const double standardDeviationCN2 = 0.1;
    const std::vector<double> priorFrequencies{ 0.1, 0.2, 0.7 };
    const std::vector<double> meanDepthValues{ 0, 1, 2 };
    const std::vector<double> likelihoodOfAllCN{ 0.1, 0.3, 0.1 };

    {
        const CopyNumberGenotyper genotyper(
            maxCopyNumber, depthScaleFactor, standardDeviationCN2, meanDepthValues, priorFrequencies);
        std::pair<int, double> bestGenotypeAndPosterior = genotyper.getBestGenotypeAndPosterior(likelihoodOfAllCN);
        EXPECT_EQ(1, bestGenotypeAndPosterior.first);
        EXPECT_EQ(0.6, bestGenotypeAndPosterior.second);
    }
}

TEST(CopyNumberGenotyping, getLikelihoodAndPvalue)
{
    const int maxCopyNumber = 2;
    const double depthScaleFactor = 1;
    const double standardDeviationCN2 = 0.1;
    const std::vector<double> priorFrequencies{ 0.1, 0.2, 0.7 };
    const std::vector<double> meanDepthValues{ 0, 1, 2 };

    {
        const CopyNumberGenotyper genotyper(
            maxCopyNumber, depthScaleFactor, standardDeviationCN2, meanDepthValues, priorFrequencies);

        {
            std::pair<double, double> likelihoodAndPvalue = genotyper.genotypeLikelihoodAndPvalue(1, 1.35);
            EXPECT_NEAR(5.39943e-06, likelihoodAndPvalue.first, 1e-10);
            EXPECT_NEAR(3.715494e-07, likelihoodAndPvalue.second, 1e-10);
        }

        {
            std::pair<double, double> likelihoodAndPvalue = genotyper.genotypeLikelihoodAndPvalue(1, 1.05);
            EXPECT_NEAR(0.8787826, likelihoodAndPvalue.first, 1e-4);
            EXPECT_NEAR(0.23975, likelihoodAndPvalue.second, 1e-4);
        }
    }
}
