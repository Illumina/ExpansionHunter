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

#include <common/ref_genome.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "thirdparty/json/json.hpp"
#include <boost/optional.hpp>

#include "common/common.h"
#include "common/genomic_region.h"
#include "region_spec/RegionBlueprint.hh"

class RegionSpec
{
public:
    RegionSpec(
        const std::string& regionId, const RegionBlueprint& regionBlueprint, AlleleCount expectedAlleleCount,
        const Region& referenceRegion);

    bool operator==(const RegionSpec& other) const;

    const std::string& regionId() const { return regionId_; }
    const Region& referenceRegion() const { return referenceRegion_; }

    const std::vector<Region>& offtargetRegions() const { return offtargetRegions_; }
    void setOfftargetRegions(const std::vector<Region>& offtargetRegions) { offtargetRegions_ = offtargetRegions; }

    const RegionBlueprint& regionBlueprint() const { return regionBlueprint_; }

    AlleleCount expectedAlleleCount() const { return expectedAlleleCount_; }

private:
    std::string regionId_;
    RegionBlueprint regionBlueprint_;
    std::vector<Region> offtargetRegions_;
    AlleleCount expectedAlleleCount_;
    Region referenceRegion_;
};

using RegionCatalog = std::map<std::string, RegionSpec>;

RegionCatalog loadRegionSpecsFromDisk(const std::string& specsPath, const RefGenome& referencePath, Sex sampleSex);
