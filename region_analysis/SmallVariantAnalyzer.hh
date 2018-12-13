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
        boost::optional<graphtools::NodeId> optionalRefNode, double haplotypeDepth)
        : VariantAnalyzer(std::move(variantId), expectedAlleleCount, graph, std::move(nodeIds))
        , variantSubtype_(variantSubtype)
        , haplotypeDepth_(haplotypeDepth)
        , optionalRefNode_(optionalRefNode)
        , alignmentClassifier_(nodeIds_)
        , smallVariantGenotyper_(haplotypeDepth_, expectedAlleleCount_)
        , allelePresenceChecker_(haplotypeDepth_)
    {
        verboseLogger_ = spdlog::get("verbose");
        // Only indels are allowed
        assert(nodeIds_.size() <= 2);
    }

    ~SmallVariantAnalyzer() = default;

    std::unique_ptr<VariantFindings> analyze() const override;

    void processMates(
        const reads::Read& read, const graphtools::GraphAlignment& readAlignment, const reads::Read& mate,
        const graphtools::GraphAlignment& mateAlignment) override;

    double haplotypeDepth() const { return haplotypeDepth_; }

protected:
    int countReadsSupportingNode(graphtools::NodeId nodeId) const;

    VariantSubtype variantSubtype_;
    double haplotypeDepth_;
    boost::optional<graphtools::NodeId> optionalRefNode_;

    ClassifierOfAlignmentsToVariant alignmentClassifier_;
    SmallVariantGenotyper smallVariantGenotyper_;
    AllelePresenceChecker allelePresenceChecker_;

    std::shared_ptr<spdlog::logger> verboseLogger_;
};

}
