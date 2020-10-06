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

#pragma once

#include <memory>

#include "common/Common.hh"
#include "common/Parameters.hh"
#include "locus_spec/CnvLocusSpec.hh"
#include "sample_analysis/DepthNormalization.hh"
#include "stats/LocusStats.hh"
#include "workflow/FeatureAnalyzer.hh"
#include "workflow/ReadCounter.hh"
#include "workflow/VariantFindings.hh"

namespace ehunter
{

class CnvVariantAnalyzer : public FeatureAnalyzer
{
public:
    CnvVariantAnalyzer(
        std::string variantId, double regionLength, CnvVariantType variantType, CopyNumberBySex contigCopyNumber,
        CnvGenotyperParameters cnvParameters, std::shared_ptr<ReadCounter> feature);
    ~CnvVariantAnalyzer() override = default;

    std::vector<std::shared_ptr<Feature>> features() override;
    CnvVariantFindings analyze(DepthNormalizer genomeDepthNormalizer) const;
    const std::string& variantId() const { return variantId_; }
    CnvVariantType variantType() const { return variantType_; }
    CopyNumberBySex contigCopyNumber() const { return contigCopyNumber_; }

protected:
    std::string variantId_;
    double regionLength_;
    CnvVariantType variantType_;
    CopyNumberBySex contigCopyNumber_;

private:
    CnvGenotyperParameters cnvParameters_;
    std::shared_ptr<ReadCounter> counter_;
};
}
