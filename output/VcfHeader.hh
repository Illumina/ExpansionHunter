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

#include <iostream>
#include <map>
#include <string>

#include "region_analysis/VariantFindings.hh"
#include "region_spec/LocusSpecification.hh"

namespace ehunter
{

enum class FieldType
{
    kInfo = 0,
    kFilter = 1,
    kFormat = 2,
    kAlt = 3
};

struct FieldDescription
{
    FieldDescription(
        FieldType fieldType, std::string id, std::string number, std::string contentType, std::string description);
    FieldType fieldType;
    std::string id;
    std::string number;
    std::string contentType;
    std::string description;
};

using FieldDescriptionIdentifier = std::pair<FieldType, std::string>;
using FieldDescriptionCatalog = std::map<FieldDescriptionIdentifier, FieldDescription>;

// Generates VCF field descriptions required for a given variant call
class FieldDescriptionWriter : public VariantFindingsVisitor
{
public:
    FieldDescriptionWriter(const LocusSpecification& regionSpec, const VariantSpecification& variantSpec)
        : regionSpec_(regionSpec)
        , variantSpec_(variantSpec)
    {
    }

    ~FieldDescriptionWriter() = default;

    void visit(const RepeatFindings* repeatFindingsPtr) override;
    void visit(const SmallVariantFindings* smallVariantFindingsPtr) override;

    void tryAddingFieldDescription(
        FieldType fieldType, const std::string& id, const std::string& number, const std::string& contentType,
        const std::string& description);

    void dumpTo(FieldDescriptionCatalog& descriptionCatalog);

private:
    const LocusSpecification& regionSpec_;
    const VariantSpecification& variantSpec_;
    FieldDescriptionCatalog fieldDescriptions_;
};

void outputVcfHeader(const RegionCatalog& regionCatalog, const SampleFindings& sampleFindings, std::ostream& out);

std::ostream& operator<<(std::ostream& out, FieldType fieldType);
std::ostream& operator<<(std::ostream& out, const FieldDescription& fieldDescription);

}