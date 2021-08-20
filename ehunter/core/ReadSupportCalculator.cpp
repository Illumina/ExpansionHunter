//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "core/ReadSupportCalculator.hh"

#include <algorithm>

namespace ehunter
{

namespace
{
int countConsistentReads(const CountTable& table, int alleleSize)
{
    int readCount = 0;

    for (const auto& sizeAndCount : table)
    {
        const int size = sizeAndCount.first;
        const int count = sizeAndCount.second;

        if (size <= alleleSize)
        {
            readCount += count;
        }
    }

    return readCount;
}
}

int ReadSupportCalculator::getCountOfConsistentSpanningReads(int alleleSize) const
{
    return spanningReadCounts_.countOf(alleleSize);
}

int ReadSupportCalculator::getCountOfConsistentFlankingReads(int alleleSize) const
{
    return countConsistentReads(flankingReadCounts_, alleleSize);
}

int ReadSupportCalculator::getCountOfConsistentRepeatReads(int alleleSize) const
{
    return countConsistentReads(inrepeatReadCounts_, alleleSize);
}

}
