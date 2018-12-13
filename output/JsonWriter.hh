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

#include "region_analysis/VariantFindings.hh"
#include "region_spec/LocusSpecification.hh"

#include "thirdparty/json/json.hpp"

namespace ehunter
{

class VariantJsonWriter : public VariantFindingsVisitor
{
public:
    VariantJsonWriter(const LocusSpecification& regionSpec, const VariantSpecification& variantSpec, int readLength)
        : regionSpec_(regionSpec)
        , variantSpec_(variantSpec)
        , readLength_(readLength)
    {
    }

    ~VariantJsonWriter() = default;
    void visit(const RepeatFindings* repeatFindingsPtr);
    void visit(const SmallVariantFindings* smallVariantFindingsPtr);
    nlohmann::json record() const { return record_; }

private:
    const LocusSpecification& regionSpec_;
    const VariantSpecification& variantSpec_;
    const int readLength_;
    nlohmann::json record_;
};

class JsonWriter
{
public:
    JsonWriter(
        const std::string& sampleName, int readLength, const RegionCatalog& regionCatalog,
        const SampleFindings& sampleFindings);

    void write(std::ostream& out);

private:
    const std::string sampleName_;
    const int readLength_;
    const RegionCatalog& regionCatalog_;
    const SampleFindings& sampleFindings_;
};

std::ostream& operator<<(std::ostream& out, JsonWriter& jsonWriter);

}
