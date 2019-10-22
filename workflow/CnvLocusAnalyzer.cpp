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

#include "workflow/CnvLocusAnalyzer.hh"

#include <string>

#include "spdlog/spdlog.h"

#include "GraphVariantAnalyzer.hh"
#include "genotyping/CopyNumberCaller.hh"
#include "region_spec/VariantSpecification.hh"
#include "workflow/CnvVariantAnalyzer.hh"
#include "workflow/FeatureAnalyzer.hh"
#include "workflow/ReadCountAnalyzer.hh"

using std::shared_ptr;
using std::string;
using std::vector;

namespace ehunter
{

using std::static_pointer_cast;

CnvLocusAnalyzer::CnvLocusAnalyzer(double minLocusCoverage, string locusId, CnvLocusSubtype locusSubtype)
    : minLocusCoverage_(minLocusCoverage)
    , locusId_(std::move(locusId))
    , locusSubtype_(std::move(locusSubtype))
{
}

void CnvLocusAnalyzer::setStats(std::shared_ptr<ReadCountAnalyzer> statsAnalyzer)
{
    readCountAnalyzer_ = std::move(statsAnalyzer);
}

void CnvLocusAnalyzer::addAnalyzer(std::shared_ptr<CnvVariantAnalyzer> variantAnalyzer)
{
    variantAnalyzers_.push_back(std::move(variantAnalyzer));
}

LocusFindings CnvLocusAnalyzer::analyze(Sex sampleSex) const
{
    LocusFindings locusFindings;

    locusFindings.optionalStats = readCountAnalyzer_->estimate(sampleSex);

    boost::optional<int> targetCopyNumber;
    std::vector<boost::optional<int>> baselineCopyNumbers;

    for (auto& analyzerPtr : variantAnalyzers_)
    {
        CnvVariantFindings varFinding = analyzerPtr->analyze();
        const VariantSubtype variantSubtype = analyzerPtr->variantSubtype();
        if (variantSubtype == VariantSubtype::kBaseline)
        {
            baselineCopyNumbers.push_back(varFinding.copyNumberCall());
        }
        else if (variantSubtype == VariantSubtype::kTarget)
        {
            targetCopyNumber = varFinding.copyNumberCall();
        }
    }

    int expectedCopyNumber = static_cast<int>(locusFindings.optionalStats->alleleCount());
    boost::optional<int> cnvLocusCopyNumberCall;
    if (locusSubtype_ == CnvLocusSubtype::kOverlapping)
    {
        cnvLocusCopyNumberCall
            = callCopyNumberForOverlappingCnv(targetCopyNumber, baselineCopyNumbers, expectedCopyNumber);
    }
    else if (locusSubtype_ == CnvLocusSubtype::kNonoverlapping)
    {
        cnvLocusCopyNumberCall
            = callCopyNumberForNonOverlappingCnv(targetCopyNumber, baselineCopyNumbers, expectedCopyNumber);
    }
    std::unique_ptr<VariantFindings> cnvLocusFindingPtr(new CnvVariantFindings(cnvLocusCopyNumberCall));
    locusFindings.findingsForEachVariant.emplace(locusId_, std::move(cnvLocusFindingPtr));

    return locusFindings;
}

vector<shared_ptr<FeatureAnalyzer>> CnvLocusAnalyzer::featureAnalyzers()
{
    vector<shared_ptr<FeatureAnalyzer>> features;
    for (const auto& variant : variantAnalyzers_)
    {
        features.push_back(variant);
    }

    if (readCountAnalyzer_ != nullptr)
    {
        features.push_back(static_pointer_cast<FeatureAnalyzer>(readCountAnalyzer_));
    }

    return features;
}
}