//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"
#include "thirdparty/json/json.hpp"

#include "common/Common.hh"
#include "common/GenomicRegion.hh"
#include "common/Reference.hh"
#include "region_spec/VariantSpecification.hh"

namespace ehunter
{

using RegionId = std::string;

class LocusSpecification
{
public:
    LocusSpecification(
        RegionId regionId, std::vector<Region> referenceLoci, AlleleCount expectedAlleleCount,
        graphtools::Graph regionGraph);

    const RegionId& regionId() const { return regionId_; }
    /*
     * List of all regions in the reference this graph describes
     * i.e. where we expect relevant reads to align
     */
    const std::vector<Region>& referenceLoci() const { return referenceLoci_; }
    /*
     * List of regions that additional relevant reads might be found
     * Require filtering or special considerations
     */
    const std::vector<Region>& offtargetLoci() const { return offtargetLoci_; }
    void setOfftargetLoci(const std::vector<Region>& offtargetLoci) { offtargetLoci_ = offtargetLoci; }
    const graphtools::Graph& regionGraph() const { return regionGraph_; }
    AlleleCount expectedAlleleCount() const { return expectedAlleleCount_; }

    const std::vector<VariantSpecification>& variantSpecs() const { return variantSpecs_; }
    void addVariantSpecification(
        std::string id, VariantClassification classification, Region referenceLocus,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode);

    const VariantSpecification& getVariantSpecById(const std::string& variantSpecId) const;

private:
    std::string regionId_;
    std::vector<Region> referenceLoci_;
    AlleleCount expectedAlleleCount_;
    std::vector<Region> offtargetLoci_;
    graphtools::Graph regionGraph_;
    std::vector<VariantSpecification> variantSpecs_;
};

using RegionCatalog = std::map<RegionId, LocusSpecification>;

}
