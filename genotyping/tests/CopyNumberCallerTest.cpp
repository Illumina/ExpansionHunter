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

#include "genotyping/CopyNumberCaller.hh"
#include "gtest/gtest.h"
#include <numeric>
#include <vector>

#include "common/Common.hh"

using boost::optional;
using namespace ehunter;

TEST(CopyNumberCalling, NonOverlappingCNV)
{
    std::vector<boost::optional<int>> baselineCopyNumbers{ 2 };
    // For non-overlapping CNVs
    // Target CN is no-call
    EXPECT_EQ(boost::none, callCopyNumberForNonOverlappingCNV(boost::none, baselineCopyNumbers, 2));
    // Baseline CN equals expected CN
    EXPECT_EQ(1, *callCopyNumberForNonOverlappingCNV(3, baselineCopyNumbers, 2));
    EXPECT_EQ(-2, *callCopyNumberForNonOverlappingCNV(0, baselineCopyNumbers, 2));
    baselineCopyNumbers = { 2, 2 };
    EXPECT_EQ(1, *callCopyNumberForNonOverlappingCNV(3, baselineCopyNumbers, 2));
    EXPECT_EQ(-2, *callCopyNumberForNonOverlappingCNV(0, baselineCopyNumbers, 2));
    // Baseline CN has no-call
    baselineCopyNumbers = { 2, boost::none };
    EXPECT_EQ(1, *callCopyNumberForNonOverlappingCNV(3, baselineCopyNumbers, 2));
    EXPECT_EQ(-2, *callCopyNumberForNonOverlappingCNV(0, baselineCopyNumbers, 2));
    // Baseline CNs don't agree
    baselineCopyNumbers = { 2, 3, boost::none };
    EXPECT_EQ(boost::none, callCopyNumberForNonOverlappingCNV(3, baselineCopyNumbers, 2));
    EXPECT_EQ(boost::none, callCopyNumberForNonOverlappingCNV(0, baselineCopyNumbers, 2));
    // Baseline CN equals Target CN
    baselineCopyNumbers = { 3, 3 };
    EXPECT_EQ(0, *callCopyNumberForNonOverlappingCNV(3, baselineCopyNumbers, 2));
    // Baseline CN is not equal to expected CN
    EXPECT_EQ(boost::none, callCopyNumberForNonOverlappingCNV(2, baselineCopyNumbers, 2));
    // Baseline CN is no-call
    baselineCopyNumbers = { boost::none, boost::none };
    EXPECT_EQ(-2, *callCopyNumberForNonOverlappingCNV(0, baselineCopyNumbers, 2));
}

TEST(CopyNumberCalling, OverlappingCNV)
{
    // For overlapping CNVs
    // Target CN is no-call
    std::vector<boost::optional<int>> baselineCopyNumbers{ 2, 2 };
    EXPECT_EQ(boost::none, callCopyNumberForOverlappingCNV(boost::none, baselineCopyNumbers, 2));
    EXPECT_EQ(1, *callCopyNumberForOverlappingCNV(3, baselineCopyNumbers, 2));
    // Baseline CN is not equal to expected CN. But for overlapping CNV it does not matter
    baselineCopyNumbers = { 3, 3 };
    EXPECT_EQ(-1, *callCopyNumberForOverlappingCNV(2, baselineCopyNumbers, 2));
    // Absolute CN has to be >= 0
    EXPECT_EQ(boost::none, callCopyNumberForOverlappingCNV(0, baselineCopyNumbers, 2));
    // Baseline CN has no-call
    baselineCopyNumbers = { boost::none, 2 };
    EXPECT_EQ(boost::none, callCopyNumberForOverlappingCNV(3, baselineCopyNumbers, 2));
    // Baseline CNs don't agree
    baselineCopyNumbers = { 2, 3 };
    EXPECT_EQ(boost::none, callCopyNumberForOverlappingCNV(3, baselineCopyNumbers, 2));
    // Baseline CN is no-call
    baselineCopyNumbers = { boost::none, boost::none };
    EXPECT_EQ(boost::none, callCopyNumberForOverlappingCNV(2, baselineCopyNumbers, 2));
}