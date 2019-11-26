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

#include "workflow/SmnLocusAnalyzer.hh"

#include <string>

#include "spdlog/spdlog.h"

#include "genotyping/CopyNumberCaller.hh"
#include "workflow/CnvVariantAnalyzer.hh"
#include "workflow/FeatureAnalyzer.hh"
#include "workflow/LinearSmallVariantAnalyzer.hh"
#include "workflow/ReadCountAnalyzer.hh"

using std::shared_ptr;
using std::string;
using std::vector;

namespace ehunter
{

using std::static_pointer_cast;

SmnLocusAnalyzer::SmnLocusAnalyzer(string locusId, std::vector<ParalogOutputVariant> outputVariants)
    : ParalogLocusAnalyzer(locusId, outputVariants)
{
}

static int findMode(std::vector<boost::optional<int>> calls)
{
    std::vector<int> histogram(10, 0);
    for (int i = 0; i != 10; ++i)
    {
        if (calls[i])
        {
            ++histogram[*calls[i]];
        }
    }
    return std::max_element(histogram.begin(), histogram.end()) - histogram.begin();
}

LocusFindings SmnLocusAnalyzer::analyze(Sex sampleSex, boost::optional<DepthNormalizer> genomeDepthNormalizer)
{
    LocusFindings locusFindings;

    locusFindings.optionalStats = readCountAnalyzer_->estimate(sampleSex);

    updateVariantFindings(genomeDepthNormalizer);
    std::vector<boost::optional<int>> copyNumberCalls;
    for (auto finding : smallVariantFindings_)
    {
        auto copyNumberPerBase = finding.copyNumber();
        if (copyNumberPerBase)
        {
            auto callPerBase = *copyNumberPerBase;
            copyNumberCalls.push_back(callPerBase.first);
        }
    }
    int smn1CopyNumberCall = findMode(copyNumberCalls);
    auto outputVariantId = "SMN1";
    std::unique_ptr<VariantFindings> cnvLocusFindingPtr(
        new CnvVariantFindings(outputVariantId, smn1CopyNumberCall, smn1CopyNumberCall));
    locusFindings.findingsForEachVariant.emplace(locusId_, std::move(cnvLocusFindingPtr));

    return locusFindings;
}
}