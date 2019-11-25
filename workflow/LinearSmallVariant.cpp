//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Xiao Chen <xchen2@illumina.com>,
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

#include "workflow/LinearSmallVariant.hh"

#include "workflow/LinearModel.hh"

using std::shared_ptr;
using std::vector;

namespace ehunter
{

LinearSmallVariant::LinearSmallVariant(
    std::shared_ptr<LinearModel> model, SmallVariantLocations locations, SmallVariantBases bases, boost::optional<int> mapqCutoff)
    : model_(std::move(model))
    , locations_(std::move(locations))
    , bases_(std::move(bases))
    , mapqCutoff_(std::move(mapqCutoff))
{
}

shared_ptr<RegionModel> LinearSmallVariant::model() { return model_; }

void LinearSmallVariant::summarize(const MappedRead& read)
{
    if (mapqCutoff_)
    {
        if (read.mapq() >= *mapqCutoff_)
        {
            ++numGeneAReads_;
            ++numGeneBReads_;
        }
    }
}

void LinearSmallVariant::summarize(const MappedRead& read, const MappedRead& mate)
{
    summarize(read);
    summarize(mate);
}
}