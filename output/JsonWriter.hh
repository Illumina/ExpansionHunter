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

#include "common/Parameters.hh"
#include "region_analysis/VariantFindings.hh"
#include "region_spec/LocusSpecification.hh"

#include "thirdparty/json/json.hpp"

namespace ehunter
{

class VariantJsonWriter : public VariantFindingsVisitor
{
public:
    VariantJsonWriter(
        const SampleParameters& sampleParams,
        const ReferenceContigInfo& contigInfo,
        const LocusSpecification& regionSpec,
        const VariantSpecification& variantSpec)
        : sampleParams_(sampleParams)
        , contigInfo_(contigInfo)
        , regionSpec_(regionSpec)
        , variantSpec_(variantSpec)
    {
    }

    ~VariantJsonWriter() = default;
    void visit(const RepeatFindings* repeatFindingsPtr);
    void visit(const SmallVariantFindings* smallVariantFindingsPtr);
    nlohmann::json record() const { return record_; }

private:
    const SampleParameters& sampleParams_;
    const ReferenceContigInfo& contigInfo_;
    const LocusSpecification& regionSpec_;
    const VariantSpecification& variantSpec_;
    nlohmann::json record_;
};

class JsonWriter
{
public:
    JsonWriter(
        const SampleParameters& sampleParams,
        const ReferenceContigInfo& contigInfo,
        const RegionCatalog& regionCatalog,
        const SampleFindings& sampleFindings);

    void write(std::ostream& out);

private:
    const SampleParameters& sampleParams_;
    const ReferenceContigInfo& contigInfo_;
    const RegionCatalog& regionCatalog_;
    const SampleFindings& sampleFindings_;
};

std::ostream& operator<<(std::ostream& out, JsonWriter& jsonWriter);

}
