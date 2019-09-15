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

using std::vector;

namespace ehunter
{

bool isFullyContained(int readContig, int64_t readStart, int64_t readEnd, const GenomicRegion& region)
{
    return ((readContig == region.contigIndex()) && (region.start() <= readStart && readEnd <= region.end()));
}

void dispatch(const MappedRead& read, const MappedRead& mate, const std::unordered_set<RegionModel*>& models)
{
    for (auto model : models)
    {
        bool readIsFullyContained = false;
        bool mateIsFullyContained = false;
        for (const auto& region : model->readExtractionRegions())
        {
            if (isFullyContained(read.contigIndex(), read.pos(), read.approximateEnd(), region))
            {
                readIsFullyContained = true;
            }

            if (isFullyContained(mate.contigIndex(), mate.pos(), mate.approximateEnd(), region))
            {
                mateIsFullyContained = true;
            }
        }

        if (readIsFullyContained && mateIsFullyContained)
        {
            // model->processFullHit(read, mate);
            model->analyze(read, mate);
        }
        else if (readIsFullyContained)
        {
            // model->processPartialHit(read, mate);
            model->analyze(read, boost::none);
        }
        else if (mateIsFullyContained)
        {
            // model->processPartialHit(mate, read);
            model->analyze(mate, boost::none);
        }
    }
}

}
