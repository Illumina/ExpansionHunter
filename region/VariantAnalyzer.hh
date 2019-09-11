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

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "region/GraphModel.hh"
#include "region/VariantFindings.hh"
#include "stats/LocusStats.hh"

namespace ehunter
{

class VariantAnalyzer
{
public:
    using SPtr = std::shared_ptr<VariantAnalyzer>;

    explicit VariantAnalyzer(std::string variantId);
    virtual ~VariantAnalyzer() = default;

    const std::string& variantId() const { return variantId_; }
    virtual std::unique_ptr<VariantFindings> analyze(const LocusStats& stats) const = 0;

    const std::vector<GraphFeature::SPtr>& regionFeaturePtrs() const { return regionFeaturePtrs_; }

private:
    std::string variantId_;
    // TODO: Introduce general "features"
    std::vector<GraphFeature::SPtr> regionFeaturePtrs_;
};

}
