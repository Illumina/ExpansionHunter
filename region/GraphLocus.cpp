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

#include "region/GraphLocus.hh"

#include <string>

using std::string;

namespace ehunter
{

LocusFindings GraphLocus::analyze(Sex /*sampleSex*/) const
{
    LocusFindings locusFindings;

    // locusFindings.optionalStats = statsCalculator_.estimate(sampleSex);

    // if (locusFindings.optionalStats
    //    && locusFindings.optionalStats->depth() >= locusSpec().genotyperParameters().minLocusCoverage)
    //{
    for (auto& variantAnalyzerPtr : variantAnalyzerPtrs_)
    {
        const LocusStats& locusStats = *locusFindings.optionalStats;
        std::unique_ptr<VariantFindings> variantFindingsPtr = variantAnalyzerPtr->analyze(locusStats);
        const string& variantId = variantAnalyzerPtr->variantId();
        locusFindings.findingsForEachVariant.emplace(variantId, std::move(variantFindingsPtr));
    }
    //}

    return locusFindings;
}

}