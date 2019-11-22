//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

#include "sample_analysis/ReadDispatch.hh"

using std::unordered_set;
using std::vector;

namespace ehunter
{

bool isFullyContained(const MappedRead& read, const GenomicRegion& region)
{
    return (
        (read.contigIndex() == region.contigIndex())
        && (region.start() <= read.pos() && read.approximateEnd() <= region.end()));
}

bool checkIfMapNearby(const MappedRead& read, const MappedRead& mate)
{
    const int kMaxDistance = 1000;
    return (read.contigIndex() == mate.contigIndex()) && (std::abs(read.pos() - mate.pos()) < kMaxDistance);
}

void dispatch(const MappedRead& read, const MappedRead& mate, const unordered_set<RegionModel*>& models)
{
    for (auto model : models)
    {
        bool readIsFullyContained = false;
        bool mateIsFullyContained = false;
        for (const auto& region : model->readExtractionRegions())
        {
            if (isFullyContained(read, region))
            {
                readIsFullyContained = true;
            }

            if (isFullyContained(mate, region))
            {
                mateIsFullyContained = true;
            }
        }

        if (readIsFullyContained && mateIsFullyContained)
        {
            model->analyze(read, mate);
        }
        else if (!checkIfMapNearby(read, mate) && (readIsFullyContained || mateIsFullyContained))
        {
            model->analyze(read, mate);
        }
        else if (readIsFullyContained)
        {
            model->analyze(read);
        }
        else if (mateIsFullyContained)
        {
            model->analyze(mate);
        }
    }
}

void dispatch(const MappedRead& read, const unordered_set<RegionModel*>& models)
{
    for (auto model : models)
    {
        for (const auto& region : model->readExtractionRegions())
        {
            if (isFullyContained(read, region))
            {
                model->analyze(read);
            }
        }
    }
}

}
