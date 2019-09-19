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
#include <math.h>
#include <numeric>
#include <set>

using boost::optional;
using std::set;
using std::vector;

namespace ehunter
{
boost::optional<int> callCopyNumber(
    const std::vector<boost::optional<int>>& baselineCopyNumbers, boost::optional<int> targetCopyNumber,
    bool baselineExpectedNormal, int expectedBaselineCopyNumber)
{
    if (targetCopyNumber == boost::none)
    {
        return boost::none;
    }

    const std::set<boost::optional<int>> baselineCopyNumberSet(baselineCopyNumbers.begin(), baselineCopyNumbers.end());
    const int baselineCopyNumberSetSize = (int)baselineCopyNumberSet.size();
    const bool hasNoCall = baselineCopyNumberSet.find(boost::none) != baselineCopyNumberSet.end();
    // For overlapping CNVs. Baseline no-calls are not allowed. No need to match expected baseline CN
    if (!baselineExpectedNormal)
    {
        if (baselineCopyNumberSetSize == 1 && !hasNoCall)
        {
            auto firstElement = baselineCopyNumberSet.begin();
            return *targetCopyNumber - **firstElement;
        }

        return boost::none;
    }
    // For non-overlapping CNVs. Baseline no-calls are allowed.
    else
    {
        // Baseline CN is no-call, use expected baseline CN
        if (baselineCopyNumberSetSize == 1 && hasNoCall)
        {
            return *targetCopyNumber - expectedBaselineCopyNumber;
        }
        // Baseline CNs have to agree with each other
        // and match either expected baseline CN or target CN
        if ((baselineCopyNumberSetSize == 1) || (baselineCopyNumberSetSize == 2 && hasNoCall))
        {
            auto nonNoCallElement = baselineCopyNumberSet.begin();
            if (*nonNoCallElement == boost::none)
            {
                std::advance(nonNoCallElement, 1);
            }
            const int baselineCopyNumber = **nonNoCallElement;
            if (baselineCopyNumber == expectedBaselineCopyNumber || baselineCopyNumber == *targetCopyNumber)
            {
                return *targetCopyNumber - baselineCopyNumber;
            }

            return boost::none;
        }

        return boost::none;
    }
}
}
