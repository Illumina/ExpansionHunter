//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>
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

#include "GenomeMask.hh"

#include <cassert>

namespace ehunter
{

namespace
{
inline size_t binPos(int64_t pos)
{
    static const int binSizeLog2 = 10;
    return pos >> binSizeLog2;
}
}

bool GenomeMask::query(int32_t contigId, int64_t pos) const
{
    if (contigId >= static_cast<int>(mask_.size()))
    {
        return false;
    }
    const contigMask& cmask = mask_[contigId];
    return ((binPos(pos) < cmask.size()) && cmask[binPos(pos)]);
}

void GenomeMask::addRegion(int32_t contigId, int64_t start, int64_t stop)
{
    assert(start <= stop);
    assert(start >= 0);
    assert(contigId >= 0);

    if (contigId >= static_cast<int>(mask_.size()))
    {
        mask_.resize(contigId + 1);
    }
    contigMask& cmask = mask_[contigId];
    const auto stopBin = binPos(stop);
    if (stopBin >= cmask.size())
    {
        cmask.resize(stopBin + 1);
    }
    for (auto bin = binPos(start); bin <= stopBin; ++bin)
    {
        cmask[bin] = true;
    }
}

}
