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

class ParalogLocusAnalyzer : public LocusAnalyzer
{
public:
    ParalogLocusAnalyzer(std::string locusId, ParalogOutputVariant outputVariant);
    ~ParalogLocusAnalyzer() override = default;

    const std::string& locusId() const override { return locusId_; }
    ParalogOutputVariant outputVariant() const { return outputVariant_; }
    void setStats(std::shared_ptr<ReadCountAnalyzer> statsAnalyzer);
    void addCnvAnalyzer(std::shared_ptr<CnvVariantAnalyzer> variantAnalyzer);
    LocusFindings analyze(Sex sampleSex, boost::optional<DepthNormalizer> genomeDepthNormalizer) const override;
    std::vector<std::shared_ptr<FeatureAnalyzer>> featureAnalyzers() override;

private:
    std::string locusId_;
    ParalogOutputVariant outputVariant_;
    std::shared_ptr<ReadCountAnalyzer> readCountAnalyzer_;
    std::vector<std::shared_ptr<CnvVariantAnalyzer>> variantAnalyzers_;
};
}
