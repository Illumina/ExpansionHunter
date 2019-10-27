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

#include "LocusSpecification.hh"
#include "common/Common.hh"
#include "common/GenomicRegion.hh"
#include "common/Reference.hh"
#include "locus_spec/VariantSpecification.hh"

namespace ehunter
{

using NodeToRegionAssociation = std::unordered_map<graphtools::NodeId, GenomicRegion>;

class GraphLocusSpecification : public LocusSpecification
{
public:
    GraphLocusSpecification(
        std::string locusId, CopyNumberBySex contigCopyNumber, GenomicRegion locusLocation,
        std::vector<GenomicRegion> targetReadExtractionRegions, graphtools::Graph regionGraph,
        NodeToRegionAssociation referenceRegions, GenotyperParameters genotyperParams)
        : LocusSpecification(locusId, contigCopyNumber, locusLocation, genotyperParams)
        , targetReadExtractionRegions_(std::move(targetReadExtractionRegions))
        , regionGraph_(std::move(regionGraph))
        , referenceRegions_(std::move(referenceRegions))
    {
    }
    ~GraphLocusSpecification() override = default;

    /*
     * List of all regions in the reference this graph describes
     * i.e. where we expect relevant reads to align
     */
    const std::vector<GenomicRegion>& targetReadExtractionRegions() const { return targetReadExtractionRegions_; }
    /*
     * List of regions where additional relevant reads might be found
     * Require filtering or special considerations
     */
    const std::vector<GenomicRegion>& offtargetReadExtractionRegions() const { return offtargetReadExtractionRegions_; }
    void setOfftargetReadExtractionRegions(const std::vector<GenomicRegion>& offtargetReadExtractionRegions)
    {
        offtargetReadExtractionRegions_ = offtargetReadExtractionRegions;
    }

    const graphtools::Graph& regionGraph() const { return regionGraph_; }
    const std::vector<VariantSpecification>& variantSpecs() const { return variantSpecs_; }

    void addVariantSpecification(
        std::string id, VariantClassification classification, GenomicRegion referenceLocus,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode);

    const NodeToRegionAssociation& referenceProjectionOfNodes() const { return referenceRegions_; }

private:
    std::vector<GenomicRegion> targetReadExtractionRegions_;
    std::vector<GenomicRegion> offtargetReadExtractionRegions_;
    graphtools::Graph regionGraph_;
    NodeToRegionAssociation referenceRegions_;
};
}
