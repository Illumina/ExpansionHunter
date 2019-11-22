//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
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

#include "workflow/ReadCounter.hh"

#include "workflow/LinearModel.hh"

using std::shared_ptr;
using std::vector;

namespace ehunter
{

ReadCounter::ReadCounter(
    shared_ptr<LinearModel> model, vector<GenomicRegion> targetRegions, boost::optional<int> mapqCutoff)
    : model_(std::move(model))
    , targetRegions_(std::move(targetRegions))
    , mapqCutoff_(std::move(mapqCutoff))
{
}

shared_ptr<RegionModel> ReadCounter::model() { return model_; }

int ReadCounter::getReadLength() const
{
    if (numReads_ == 0)
    {
        return 0;
    }

    return static_cast<int>(totalReadLength_ / numReads_);
}

double ReadCounter::getDepth() const
{
    const int readLength = getReadLength();
    std::int64_t numberOfStartPositions = 0;

    for (const auto& region : targetRegions_)
    {
        numberOfStartPositions += region.length() - readLength;
    }

    assert(numberOfStartPositions > 0);
    const double depth = readLength * (static_cast<double>(numReads()) / numberOfStartPositions);

    return depth;
}

void ReadCounter::summarize(const MappedRead& read)
{
    ++numReads_;
    totalReadLength_ += read.sequence().length();

    bool readIsInRegion = false;
    for (GenomicRegion region : targetRegions_)
    {
        if ((read.pos() < region.end()) & (read.pos() >= region.start()))
        {
            readIsInRegion = true;
        }
    }

    if (readIsInRegion)
    {
        if (mapqCutoff_)
        {
            if (read.mapq() >= *mapqCutoff_)
            {
                ++numReadsForCnvCounting_;
            }
        }
    }
}

void ReadCounter::summarize(const MappedRead& read, const MappedRead& mate)
{
    summarize(read);
    summarize(mate);
}
}