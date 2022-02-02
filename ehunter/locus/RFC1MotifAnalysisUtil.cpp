//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Chris Saunders <csaunders@illumina.com>
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

#include "RFC1MotifAnalysisUtil.hh"

#include <algorithm>

namespace ehunter
{

namespace
{

/// \brief Determine the last usable base in the read presented in cycle order
///
/// The last usable base is the last base occurring in the window which immediately proceeds the first non high-quality
/// window in the read. The window size is 10 bases, windows are considered high quality when all 10 bases are quantized
/// to the high-quality state.
///
/// \param[in] binaryQuals Quality vector for the read, reduced to 2 {low, high} quality states
///
/// \return Zero-indexed position of last usable base in cycle coordinates, or none if no usable bases are found
///
boost::optional<unsigned> findLastUsableReadCycle(const std::vector<uint8_t>& binaryQuals)
{
    // Size of the quality assessment window
    const unsigned qualWinSize(10);

    // Min number of high-quality bases in the window
    const unsigned minQsum(10);

    const unsigned readSize(binaryQuals.size());
    assert(readSize >= qualWinSize);
    const unsigned maxCycleIndex(readSize - (qualWinSize - 1));
    unsigned cycleIndex(0);
    for (; cycleIndex < maxCycleIndex; ++cycleIndex)
    {
        uint16_t qsum(0);
        for (unsigned winIndex(0); winIndex < qualWinSize; ++winIndex)
        {
            qsum += binaryQuals[cycleIndex + winIndex];
        }
        if (qsum < minQsum)
        {
            if (cycleIndex == 0)
            {
                return {};
            }
            else
            {
                break;
            }
        }
    }

    return (cycleIndex + qualWinSize - 2);
}

}

std::string getMinRotation(std::string str)
{
    auto minStr(str);
    const unsigned rotationCount(str.size() - 1);
    for (unsigned rotIndex(0); rotIndex < rotationCount; rotIndex++)
    {
        std::rotate(str.begin(), str.begin() + 1, str.end());
        if (str < minStr)
        {
            minStr = str;
        }
    }
    return minStr;
}

boost::optional<std::pair<unsigned, unsigned>>
findUsableReadBaseRange(std::vector<uint8_t> binaryQuals, const bool isReversed)
{
    if (isReversed)
    {
        std::reverse(binaryQuals.begin(), binaryQuals.end());
    }

    const auto cycleIndex = findLastUsableReadCycle(binaryQuals);
    if (not cycleIndex)
    {
        return {};
    }

    if (isReversed)
    {
        const unsigned readSize(binaryQuals.size());
        return std::make_pair(readSize - (*cycleIndex + 1), readSize - 1);
    }
    else
    {
        return std::make_pair(0u, *cycleIndex);
    }
}

}
