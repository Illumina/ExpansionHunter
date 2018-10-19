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

#include "region_analysis/RepeatFindings.hh"
#include "region_spec/RegionBlueprint.hh"
#include "region_spec/RegionSpec.hh"

#include "thirdparty/json/json.hpp"

class JsonWriter
{
public:
    JsonWriter(
        const std::string& sampleName, int readLength, const RegionCatalog& regionSpecs,
        const SampleFindings& sampleFindings);

    void write(std::ostream& out);

private:
    void addRegionFindings(const std::string& regionId, const RegionFindings& regionFindings, nlohmann::json& array);
    void addRepeatFindings(
        const RegionBlueprintComponent& repeatBlueprint, const RepeatFindings& repeatFindings, nlohmann::json& array);

    std::string encodeRepeatAlleleSupport(const std::string& repeatUnit, const RepeatFindings& repeatFindings);

    const std::string sampleName_;
    const int readLength_;
    const RegionCatalog& regionSpecs_;
    const SampleFindings& sampleFindings_;
};

std::ostream& operator<<(std::ostream& out, JsonWriter& jsonWriter);
