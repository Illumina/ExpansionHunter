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

static int findMode(std::vector<boost::optional<int>> copyNumberCalls)
{
    int maxCopyNumber = 0;
    for (auto callPerSite : copyNumberCalls)
    {
        if (callPerSite)
        {
            if (*callPerSite > maxCopyNumber)
            {
                maxCopyNumber = *callPerSite;
            }
        }
    }
    std::vector<int> callsCount(maxCopyNumber, 0);
    for (int i = 0; i != maxCopyNumber; ++i)
    {
        if (copyNumberCalls[i])
        {
            ++callsCount[*copyNumberCalls[i]];
        }
    }
    return std::max_element(callsCount.begin(), callsCount.end()) - callsCount.begin();
}

LocusFindings SmnLocusAnalyzer::analyze(Sex sampleSex, boost::optional<DepthNormalizer> genomeDepthNormalizer)
{
    LocusFindings locusFindings;

    locusFindings.optionalStats = readCountAnalyzer_->estimate(sampleSex);

    updateVariantFindings(genomeDepthNormalizer);
    auto outputVariantId1 = "SMN1";
    auto outputVariantId2 = "SMN2";

    boost::optional<int> totalCopyNumber;
    boost::optional<int> intactCopyNumber;
    boost::optional<int> smn1CopyNumberCall;
    boost::optional<int> smn2CopyNumberCall;
    for (auto finding : cnvFindings_)
    {
        if (finding.variantId() == "Exon1-6")
        {
            totalCopyNumber = finding.absoluteCopyNumber();
        }
        if (finding.variantId() == "Exon7-8")
        {
            intactCopyNumber = finding.absoluteCopyNumber();
        }
    }
    /*
    if (totalCopyNumber && intactCopyNumber)
    {
        int SmnDelta = *totalCopyNumber - *intactCopyNumber;
    }
    */

    if (intactCopyNumber)
    {
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
        smn1CopyNumberCall = findMode(copyNumberCalls);
        smn2CopyNumberCall = *intactCopyNumber - *smn1CopyNumberCall;
    }

    std::unique_ptr<VariantFindings> Smn1FindingPtr(
        new CnvVariantFindings(outputVariantId1, smn1CopyNumberCall, smn1CopyNumberCall));
    locusFindings.findingsForEachVariant.emplace(outputVariantId1, std::move(Smn1FindingPtr));
    std::unique_ptr<VariantFindings> Smn2FindingPtr(
        new CnvVariantFindings(outputVariantId2, smn2CopyNumberCall, smn2CopyNumberCall));
    locusFindings.findingsForEachVariant.emplace(outputVariantId2, std::move(Smn2FindingPtr));

    return locusFindings;
}
}