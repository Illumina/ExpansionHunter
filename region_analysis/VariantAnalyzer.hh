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

#include <memory>
#include <string>
#include <vector>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "common/Common.hh"
#include "reads/Read.hh"
#include "region_analysis/VariantFindings.hh"

namespace ehunter
{

class VariantAnalyzer
{
public:
    VariantAnalyzer(
        std::string variantId, AlleleCount expectedAlleleCount, const graphtools::Graph& graph,
        std::vector<graphtools::NodeId> nodeIds)
        : variantId_(std::move(variantId))
        , expectedAlleleCount_(expectedAlleleCount)
        , graph_(graph)
        , nodeIds_(std::move(nodeIds))
    {
    }
    virtual ~VariantAnalyzer() = default;

    virtual void processMates(
        const reads::Read& read, const graphtools::GraphAlignment& readAlignment, const reads::Read& mate,
        const graphtools::GraphAlignment& mateAlignment)
        = 0;

    virtual std::unique_ptr<VariantFindings> analyze() const = 0;

    const std::string& variantId() const { return variantId_; }
    AlleleCount expectedAlleleCount() const { return expectedAlleleCount_; }
    const graphtools::Graph& graph() const { return graph_; }
    const std::vector<graphtools::NodeId>& nodeIds() const { return nodeIds_; }

protected:
    std::string variantId_;
    AlleleCount expectedAlleleCount_;
    const graphtools::Graph& graph_;
    std::vector<graphtools::NodeId> nodeIds_;
};

}
