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

#include "region/LocusAnalyzer.hh"

#include <unordered_set>

using std::unordered_set;
using std::vector;

namespace ehunter
{

vector<Region::SPtr> extractRegionModels(const vector<LocusAnalyzer::SPtr>& locusAnalyzerPtrs)
{
    unordered_set<Region::SPtr> regionPtrs;

    for (const auto& locusAnalyzerPtr : locusAnalyzerPtrs)
    {
        for (const auto& variantAnalyzerPtr : locusAnalyzerPtr->variantAnalyzerPtrs())
        {
            for (const auto& regionFeaturePtr : variantAnalyzerPtr->regionFeaturePtrs())
            {
                regionPtrs.insert(regionFeaturePtr->regionModelPtr());
            }
        }
    }

    return vector<Region::SPtr>(regionPtrs.begin(), regionPtrs.end());
}

}
