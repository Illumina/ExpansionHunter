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

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"

#include "common/Parameters.hh"
#include "genotyping/AlleleChecker.hh"
#include "region_spec/VariantSpecification.hh"
#include "workflow/GraphVariantAnalyzer.hh"

namespace ehunter
{

class SmallVariantFeature;

class SmallVariantAnalyzer : public GraphVariantAnalyzer
{
public:
    SmallVariantAnalyzer(
        std::shared_ptr<SmallVariantFeature> smallVariantFeature, std::string variantId, VariantSubtype variantSubtype,
        boost::optional<graphtools::NodeId> optionalRefNode);
    ~SmallVariantAnalyzer() override = default;

    std::vector<std::shared_ptr<ModelFeature>> features() override;
    std::unique_ptr<VariantFindings> analyze(const LocusStats& stats) const override;

private:
    std::shared_ptr<SmallVariantFeature> smallVariantFeature_;
    VariantSubtype variantSubtype_;
    boost::optional<graphtools::NodeId> optionalRefNode_;
    GenotyperParameters genotyperParams_;
    AlleleChecker allelePresenceChecker_;
};

}
