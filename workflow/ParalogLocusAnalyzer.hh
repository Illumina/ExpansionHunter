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

#include <memory>
#include <string>
#include <vector>

#include "locus_spec/ParalogLocusSpec.hh"
#include "workflow/LocusAnalyzer.hh"

namespace ehunter
{

class CnvVariantAnalyzer;
class ReadCountAnalyzer;
class LinearSmallVariantAnalyzer;

class ParalogLocusAnalyzer : public LocusAnalyzer
{
public:
    ParalogLocusAnalyzer(std::string locusId, std::vector<ParalogOutputVariant> outputVariants);
    // ~ParalogLocusAnalyzer() override = default;
    virtual ~ParalogLocusAnalyzer() = default;

    const std::string& locusId() const override { return locusId_; }
    std::vector<ParalogOutputVariant> outputVariants() const { return outputVariants_; }
    void setStats(std::shared_ptr<ReadCountAnalyzer> statsAnalyzer);
    void addCnvAnalyzer(std::shared_ptr<CnvVariantAnalyzer> variantAnalyzer);
    void addSmallVariantAnalyzer(std::shared_ptr<LinearSmallVariantAnalyzer> variantAnalyzer);
    void updateVariantFindings(boost::optional<DepthNormalizer> genomeDepthNormalizer);
    virtual LocusFindings analyze(Sex sampleSex, boost::optional<DepthNormalizer> genomeDepthNormalizer) = 0;
    std::vector<std::shared_ptr<FeatureAnalyzer>> featureAnalyzers() override;

protected:
    std::string locusId_;
    std::vector<ParalogOutputVariant> outputVariants_;
    std::shared_ptr<ReadCountAnalyzer> readCountAnalyzer_;
    std::vector<std::shared_ptr<CnvVariantAnalyzer>> cnvVariantAnalyzers_;
    std::vector<std::shared_ptr<LinearSmallVariantAnalyzer>> smallVariantAnalyzers_;
    std::vector<CnvVariantFindings> cnvFindings_;
    std::vector<ParalogSmallVariantFindings> smallVariantFindings_;
};
}
