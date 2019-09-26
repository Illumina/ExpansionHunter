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

#include "workflow/OfftargetFeature.hh"

#include "workflow/GraphModel.hh"

using std::shared_ptr;
using std::string;

namespace ehunter
{

OfftargetFeature::OfftargetFeature(shared_ptr<GraphModel> model, string motif)
    : model_(std::move(model))
    , motif_(std::move(motif))
    , weightedPurityCalculator_(motif_)
{
}

shared_ptr<RegionModel> OfftargetFeature::model() { return model_; }

void OfftargetFeature::process(const MappedRead& read, const MappedRead& mate)
{
    const bool isFirstReadInrepeat = weightedPurityCalculator_.score(read.sequence()) >= 0.90;
    const bool isSecondReadInrepeat = weightedPurityCalculator_.score(mate.sequence()) >= 0.90;

    if (isFirstReadInrepeat && isSecondReadInrepeat)
    {
        ++numIrrPairs_;
    }
}

}
