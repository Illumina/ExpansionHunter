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

#include "locus_spec/CnvLocusSpec.hh"
#include "workflow/LocusAnalyzer.hh"

namespace ehunter
{

class CnvVariantAnalyzer;
class ReadCountAnalyzer;

class CnvLocusAnalyzer : public LocusAnalyzer
{
public:
    CnvLocusAnalyzer(std::string locusId, CnvLocusType locusType, CnvOutputVariant outputVariant);
    ~CnvLocusAnalyzer() override = default;

    const std::string& locusId() const override { return locusId_; }
    CnvLocusType locusType() const { return locusType_; }
    CnvOutputVariant outputVariant() const { return outputVariant_; }
    void setStats(std::shared_ptr<ReadCountAnalyzer> statsAnalyzer);
    void addAnalyzer(std::shared_ptr<CnvVariantAnalyzer> variantAnalyzer);
    LocusFindings analyze(Sex sampleSex) const override;
    std::vector<std::shared_ptr<FeatureAnalyzer>> featureAnalyzers() override;

private:
    std::string locusId_;
    CnvLocusType locusType_;
    CnvOutputVariant outputVariant_;
    std::shared_ptr<ReadCountAnalyzer> readCountAnalyzer_;
    std::vector<std::shared_ptr<CnvVariantAnalyzer>> variantAnalyzers_;
};
}
