//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"
#include "thirdparty/json/json.hpp"

#include "LocusSpec.hh"
#include "common/Common.hh"
#include "common/GenomicRegion.hh"
#include "common/Reference.hh"
#include "locus_spec/VariantSpec.hh"

namespace ehunter
{

struct GraphLocusReferenceRegions
{
    // Regions in the reference where we expect relevant reads to align
    std::vector<GenomicRegion> regionsWithReads;
    // Regions where additional relevant reads might be found that require filtering or special considerations
    std::vector<GenomicRegion> offtargetRegionsWithReads;
    std::vector<GenomicRegion> statsRegions;
};

using NodeToRegionAssociation = std::unordered_map<graphtools::NodeId, GenomicRegion>;

struct ReferenceGraph
{
    ReferenceGraph(graphtools::Graph graph, NodeToRegionAssociation nodeLocations)
        : graph(std::move(graph))
        , nodeLocations(std::move(nodeLocations))
    {
    }

    graphtools::Graph graph;
    NodeToRegionAssociation nodeLocations;
};

class GraphLocusSpec : public LocusSpec
{
public:
    GraphLocusSpec(
        std::string locusId, CopyNumberBySex contigCopyNumber, GraphLocusReferenceRegions referenceRegions,
        ReferenceGraph referenceGraph, GenotyperParameters genotyperParams)
        : LocusSpec(locusId, contigCopyNumber, genotyperParams)
        , referenceGraph_(std::move(referenceGraph))
        , referenceRegions_(std::move(referenceRegions))
    {
    }
    ~GraphLocusSpec() override = default;

    const std::vector<GenomicRegion>& regionsWithReads() const { return referenceRegions_.regionsWithReads; }
    const std::vector<GenomicRegion>& offtargetRegionsWithReads() const
    {
        return referenceRegions_.offtargetRegionsWithReads;
    }
    const std::vector<GenomicRegion>& statsRegions() const { return referenceRegions_.statsRegions; }

    const graphtools::Graph& graph() const { return referenceGraph_.graph; }
    const NodeToRegionAssociation& nodeLocations() const { return referenceGraph_.nodeLocations; }

    const std::vector<VariantSpec>& variantSpecs() const { return variantSpecs_; }
    void addVariantSpecification(
        std::string id, VariantClassification classification, GenomicRegion referenceLocus,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode);

private:
    ReferenceGraph referenceGraph_;
    GraphLocusReferenceRegions referenceRegions_;
};
}
