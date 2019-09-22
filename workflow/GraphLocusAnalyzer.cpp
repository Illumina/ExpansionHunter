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

#include "workflow/GraphLocusAnalyzer.hh"

#include <string>

#include "GraphVariantAnalyzer.hh"
#include "workflow/FeatureAnalyzer.hh"
#include "workflow/StatsAnalyzer.hh"

using std::shared_ptr;
using std::string;
using std::vector;

namespace ehunter
{

using std::static_pointer_cast;

void GraphLocusAnalyzer::setStats(std::shared_ptr<StatsAnalyzer> statsAnalyzer)
{
    statsAnalyzer_ = std::move(statsAnalyzer);
}

void GraphLocusAnalyzer::addAnalyzer(std::shared_ptr<GraphVariantAnalyzer> variantAnalyzer)
{
    variantAnalyzers_.push_back(variantAnalyzer);
}

LocusFindings GraphLocusAnalyzer::analyze(Sex /*sampleSex*/) const
{
    LocusFindings locusFindings;

    // locusFindings.optionalStats = statsCalculator_.estimate(sampleSex);

    // if (locusFindings.optionalStats
    //    && locusFindings.optionalStats->depth() >= locusSpec().genotyperParameters().minLocusCoverage)
    //{

    for (auto& analyzerPtr : variantAnalyzers_)
    {
        const LocusStats& locusStats = *locusFindings.optionalStats;
        std::unique_ptr<VariantFindings> variantFindingsPtr = analyzerPtr->analyze(locusStats);
        const string& variantId = analyzerPtr->variantId();
        locusFindings.findingsForEachVariant.emplace(variantId, std::move(variantFindingsPtr));
    }

    //}

    return locusFindings;
}

vector<shared_ptr<FeatureAnalyzer>> GraphLocusAnalyzer::featureAnalyzers()
{
    vector<shared_ptr<FeatureAnalyzer>> features;
    for (const auto& variant : variantAnalyzers_)
    {
        features.push_back(variant);
    }

    if (statsAnalyzer_ != nullptr)
    {
        features.push_back(static_pointer_cast<FeatureAnalyzer>(statsAnalyzer_));
    }

    return features;
}

}