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

#include "workflow/ParalogLocusAnalyzer.hh"

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

ParalogLocusAnalyzer::ParalogLocusAnalyzer(string locusId, ParalogOutputVariant outputVariant)
    : locusId_(std::move(locusId))
    , outputVariant_(outputVariant)
{
}

void ParalogLocusAnalyzer::setStats(std::shared_ptr<ReadCountAnalyzer> statsAnalyzer)
{
    readCountAnalyzer_ = std::move(statsAnalyzer);
}

void ParalogLocusAnalyzer::addCnvAnalyzer(std::shared_ptr<CnvVariantAnalyzer> variantAnalyzer)
{
    cnvVariantAnalyzers_.push_back(std::move(variantAnalyzer));
}

void ParalogLocusAnalyzer::addSmallVariantAnalyzer(std::shared_ptr<LinearSmallVariantAnalyzer> variantAnalyzer)
{
    smallVariantAnalyzers_.push_back(std::move(variantAnalyzer));
}

LocusFindings ParalogLocusAnalyzer::analyze(Sex sampleSex, boost::optional<DepthNormalizer> genomeDepthNormalizer) const
{
    LocusFindings locusFindings;

    locusFindings.optionalStats = readCountAnalyzer_->estimate(sampleSex);

    boost::optional<int> totalCopyNumber;
    for (auto& analyzerPtr : cnvVariantAnalyzers_)
    {
        std::cout << analyzerPtr->variantId() << "\n";
        auto depthNormalizer = *genomeDepthNormalizer;
        CnvVariantFindings varFinding = analyzerPtr->analyze(depthNormalizer);
        if (analyzerPtr->variantType() == CnvVariantType::kTarget)
        {
            totalCopyNumber = *varFinding.absoluteCopyNumber();
        }
        std::cout << analyzerPtr->variantId() << " " << *varFinding.absoluteCopyNumber() << "\n";
        // locusFindings.findingsForEachVariant.emplace(locusId_, varFinding);
    }

    for (auto& analyzerPtr : smallVariantAnalyzers_)
    {
        ParalogSmallVariantFindings varFinding = analyzerPtr->analyze(totalCopyNumber);
        std::cout << analyzerPtr->variantId() << " " << varFinding.numGeneAReads() << " "
                  << varFinding.numGeneBReads() << "\n";
        // locusFindings.findingsForEachVariant.emplace(locusId_, varFinding);
    }

    return locusFindings;
}

vector<shared_ptr<FeatureAnalyzer>> ParalogLocusAnalyzer::featureAnalyzers()
{
    vector<shared_ptr<FeatureAnalyzer>> features;
    for (const auto& variant : cnvVariantAnalyzers_)
    {
        features.push_back(variant);
    }
    for (const auto& variant : smallVariantAnalyzers_)
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