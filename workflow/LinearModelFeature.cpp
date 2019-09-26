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

#include "workflow/LinearModelFeature.hh"

#include "workflow/LinearModel.hh"

using std::shared_ptr;
using std::vector;

namespace ehunter
{

LinearModelFeature::LinearModelFeature(shared_ptr<LinearModel> modelPtr, vector<GenomicRegion> targetRegions)
    : modelPtr_(std::move(modelPtr))
    , targetRegions_(std::move(targetRegions))
{
}

shared_ptr<RegionModel> LinearModelFeature::model() { return modelPtr_; }

int LinearModelFeature::getReadLength() const
{
    if (numReads_ == 0)
    {
        return 0;
    }

    return static_cast<int>(totalReadLength_ / numReads_);
}

double LinearModelFeature::getDepth() const
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

void LinearModelFeature::addReadInfo(int readLength)
{
    ++numReads_;
    totalReadLength_ += readLength;
}

}
