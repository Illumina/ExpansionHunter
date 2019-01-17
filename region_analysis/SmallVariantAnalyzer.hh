//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#include <boost/optional.hpp>

#include "thirdparty/spdlog/spdlog.h"

#include "classification/ClassifierOfAlignmentsToVariant.hh"
#include "genotyping/AllelePresenceChecker.hh"
#include "genotyping/SmallVariantGenotyper.hh"
#include "region_analysis/VariantAnalyzer.hh"
#include "region_spec/VariantSpecification.hh"

namespace ehunter
{

class SmallVariantAnalyzer : public VariantAnalyzer
{
public:
    SmallVariantAnalyzer(
        std::string variantId, VariantSubtype variantSubtype, AlleleCount expectedAlleleCount,
        const graphtools::Graph& graph, std::vector<graphtools::NodeId> nodeIds,
        boost::optional<graphtools::NodeId> optionalRefNode)
        : VariantAnalyzer(std::move(variantId), expectedAlleleCount, graph, std::move(nodeIds))
        , variantSubtype_(variantSubtype)
        , optionalRefNode_(optionalRefNode)
        , alignmentClassifier_(nodeIds_)
        , console_(spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console"))
    {
        // Only indels are allowed
        assert(nodeIds_.size() <= 2);
    }

    ~SmallVariantAnalyzer() = default;

    std::unique_ptr<VariantFindings> analyze(const SampleParameters& params) const override;

    void processMates(
        const Read& read, const graphtools::GraphAlignment& readAlignment, const Read& mate,
        const graphtools::GraphAlignment& mateAlignment) override;

protected:
    int countReadsSupportingNode(graphtools::NodeId nodeId) const;

    VariantSubtype variantSubtype_;
    boost::optional<graphtools::NodeId> optionalRefNode_;

    ClassifierOfAlignmentsToVariant alignmentClassifier_;
    AllelePresenceChecker allelePresenceChecker_;

    std::shared_ptr<spdlog::logger> console_;
};

}
