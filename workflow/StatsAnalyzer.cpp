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

#include "workflow/StatsAnalyzer.hh"

using std::shared_ptr;
using std::vector;

namespace ehunter
{

StatsAnalyzer::StatsAnalyzer(std::shared_ptr<CountingFeature> feature)
    : feature_(std::move(feature))
{
}

vector<shared_ptr<ModelFeature>> StatsAnalyzer::features() { return { feature_ }; }

LocusStats StatsAnalyzer::estimate(Sex sampleSex) const
{
    const int readLength = feature_->getReadLength();
    const double depth = feature_->getDepth();

    return { AlleleCount::kTwo, readLength, depth };
}

}
