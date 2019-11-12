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
#include "workflow/CnvVariantAnalyzer.hh"
#include "workflow/FeatureAnalyzer.hh"
#include "workflow/ReadCountAnalyzer.hh"

using std::shared_ptr;
using std::string;
using std::vector;

namespace ehunter
{

using std::static_pointer_cast;

CnvLocusAnalyzer::CnvLocusAnalyzer(string locusId, CnvLocusType locusType, CnvOutputVariant outputVariant)
    : locusId_(std::move(locusId))
    , locusType_(locusType)
    , outputVariant_(outputVariant)
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

LocusFindings CnvLocusAnalyzer::analyze(Sex sampleSex, boost::optional<DepthNormalizer> genomeDepthNormalizer) const
{
    LocusFindings locusFindings;

    locusFindings.optionalStats = readCountAnalyzer_->estimate(sampleSex);

    boost::optional<int> targetCopyNumber;
    std::vector<boost::optional<int>> baselineCopyNumbers;

    for (auto& analyzerPtr : variantAnalyzers_)
    {
        auto depthNormalizer = *genomeDepthNormalizer;
        CnvVariantFindings varFinding = analyzerPtr->analyze(depthNormalizer);
        auto variantType = analyzerPtr->variantType();
        if (variantType == CnvVariantType::kBaseline)
        {
            baselineCopyNumbers.push_back(varFinding.copyNumberCall());
        }
        else if (variantType == CnvVariantType::kTarget)
        {
            targetCopyNumber = varFinding.copyNumberCall();
        }
    }

    int expectedCopyNumber = static_cast<int>(locusFindings.optionalStats->alleleCount());
    boost::optional<int> cnvLocusCopyNumberCall;
    if (locusType_ == CnvLocusType::kOverlapping)
    {
        cnvLocusCopyNumberCall
            = callCopyNumberForOverlappingCnv(targetCopyNumber, baselineCopyNumbers, expectedCopyNumber);
    }
    else if (locusType_ == CnvLocusType::kNonoverlapping)
    {
        cnvLocusCopyNumberCall
            = callCopyNumberForNonOverlappingCnv(targetCopyNumber, baselineCopyNumbers, expectedCopyNumber);
    }

    auto outputVariantId = outputVariant_.id;
    std::unique_ptr<VariantFindings> cnvLocusFindingPtr(
        new CnvVariantFindings(outputVariantId, cnvLocusCopyNumberCall));
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