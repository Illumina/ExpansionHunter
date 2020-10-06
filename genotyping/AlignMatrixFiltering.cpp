//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "genotyping/AlignMatrixFiltering.hh"

namespace ehunter
{
namespace strgt
{

void filter(AlignMatrix& aligns)
{
    assert(aligns.numReads() % 2 == 0);
    int numStrReads = 0;
    int numStrReadsWithIndels = 0;
    int longestStrSize = 0;
    int longestStrSizeWithoutIndels = 0;
    for (int readIndex = 0; readIndex != aligns.numReads(); ++readIndex)
    {
        const auto& topAlign = aligns.getBestAlign(readIndex);

        longestStrSize = std::max(longestStrSize, topAlign.numMotifs());

        if (topAlign.type() != StrAlign::Type::kOutside)
        {
            ++numStrReads;
            if (topAlign.numIndels() > 0)
            {
                ++numStrReadsWithIndels;
            }
            else
            {
                longestStrSizeWithoutIndels = std::max(longestStrSizeWithoutIndels, topAlign.numMotifs());
            }
        }
    }

    const int maxOutliers = std::max(1, static_cast<int>(numStrReads * 0.20));
    if (numStrReadsWithIndels == 0 or maxOutliers < numStrReadsWithIndels)
    {
        return;
    }

    const double increase = (static_cast<double>(longestStrSize) - longestStrSizeWithoutIndels) / longestStrSize;

    if (increase < 0.1)
    {
        return;
    }

    for (int readIndex = aligns.numReads() - 1; readIndex != -1; --readIndex)
    {
        if (aligns.getBestAlign(readIndex).numIndels() == 0)
        {
            continue;
        }

        if (readIndex % 2 == 0)
        {
            aligns.remove(readIndex + 1);
            aligns.remove(readIndex);
        }
        else
        {
            aligns.remove(readIndex);
            aligns.remove(readIndex - 1);
            --readIndex;
        }
    }
}

CountTable countAligns(StrAlign::Type alignType, const AlignMatrix& aligns)
{
    CountTable countTable;
    for (int readIndex = 0; readIndex != aligns.numReads(); ++readIndex)
    {
        const auto& align = aligns.getBestAlign(readIndex);
        if (align.type() == alignType)
        {
            countTable.incrementCountOf(align.numMotifs());
        }
    }

    return countTable;
}

}
}