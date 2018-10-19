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

#include <memory>
#include <string>

#include <boost/optional.hpp>

#include "graphalign/GappedAligner.hh"
#include "graphalign/KmerIndex.hh"
#include "graphcore/Graph.hh"

#include "alignment/SoftclippingAligner.hh"
#include "filtering/OrientationPredictor.hh"
#include "reads/read.h"
#include "region_analysis/RepeatAnalyzer.hh"
#include "region_analysis/RepeatFindings.hh"
#include "region_spec/RegionSpec.hh"
#include "region_spec/region_graph.h"

#pragma once

class RegionAnalyzer
{
public:
    RegionAnalyzer(
        const RegionSpec& regionSpec, double haplotypeDepth, int32_t readLength, std::ostream& alignmentStream,
        const std::string& alignerName, GraphAlignmentHeuristicsParameters alignmentParams)
        : graph_(makeRegionGraph(regionSpec.regionBlueprint()))
        , regionSpec_(regionSpec)
        , haplotypeDepth_(haplotypeDepth)
        , alignmentStream_(alignmentStream)
        , alignmentParams_(alignmentParams)
        , orientationPredictor_(readLength, &graph_)
        , graphAligner_(&graph_, alignerName, alignmentParams_)
    {
        int32_t nodeId = 0;
        for (const auto& component : regionSpec.regionBlueprint())
        {
            if (component.type() == RegionBlueprintComponent::Type::kRepeat)
            {
                int32_t repeatUnitLength = component.sequence().length();
                int32_t maxNumUnitsInRead = std::ceil(readLength / static_cast<double>(repeatUnitLength));
                repeatAnalyzers_.emplace_back(
                    component.id(), graph_, regionSpec.expectedAlleleCount(), nodeId, haplotypeDepth_,
                    maxNumUnitsInRead);
            }
            ++nodeId;
        }
    }

    RegionAnalyzer(const RegionAnalyzer&) = delete;
    RegionAnalyzer& operator=(const RegionAnalyzer&) = delete;
    RegionAnalyzer(RegionAnalyzer&&) = default;
    RegionAnalyzer& operator=(RegionAnalyzer&&) = default;

    const std::string& regionId() const { return regionSpec_.regionId(); }
    const RegionSpec& regionSpec() const { return regionSpec_; }

    void processMates(reads::Read read1, reads::Read read2);
    bool checkIfPassesSequenceFilters(const std::string& sequence) const;
    bool checkIfPassesAlignmentFilters(const graphtools::GraphAlignment& alignment) const;
    bool checkIfPassesAlignmentFilters() const;

    RegionFindings genotype();

    bool operator==(const RegionAnalyzer& other) const;

private:
    boost::optional<GraphAlignment> alignRead(reads::Read& read) const;

    graphtools::Graph graph_;
    RegionSpec regionSpec_;
    double haplotypeDepth_;

    std::ostream& alignmentStream_;
    GraphAlignmentHeuristicsParameters alignmentParams_;
    OrientationPredictor orientationPredictor_;
    SoftclippingAligner graphAligner_;

    std::vector<RepeatAnalyzer> repeatAnalyzers_;
};

std::vector<std::unique_ptr<RegionAnalyzer>> initializeRegionAnalyzers(
    const RegionCatalog& RegionCatalog, double haplotypeDepth, int readLength, const std::string& alignerName,
    std::ostream& alignmentStream);