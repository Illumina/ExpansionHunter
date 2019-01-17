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
#include <set>
#include <string>

#include "common/Parameters.hh"
#include "common/Reference.hh"
#include "region_analysis/VariantFindings.hh"
#include "region_spec/LocusSpecification.hh"

namespace ehunter
{

class VariantVcfWriter : public VariantFindingsVisitor
{
public:
    VariantVcfWriter(
        const SampleParameters& sampleParams, Reference& reference, const LocusSpecification& regionSpec,
        const VariantSpecification& variantSpec, std::ostream& out)
        : sampleParams_(sampleParams)
        , reference_(reference)
        , regionSpec_(regionSpec)
        , variantSpec_(variantSpec)
        , out_(out)
    {
    }

    ~VariantVcfWriter() = default;
    void visit(const RepeatFindings* repeatFindingsPtr) override;
    void visit(const SmallVariantFindings* smallVariantFindingsPtr) override;

private:
    const SampleParameters& sampleParams_;
    Reference& reference_;
    const LocusSpecification& regionSpec_;
    const VariantSpecification& variantSpec_;
    std::ostream& out_;
};

// TODO: Document the code after multi-unit repeat format is finalized (GT-598)
class VcfWriter
{
public:
    VcfWriter(
        const SampleParameters& sampleParams, Reference& reference, const RegionCatalog& regionCatalog,
        const SampleFindings& sampleFindings);

    friend std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter);

private:
    void writeHeader(std::ostream& out);
    void writeBody(std::ostream& out);
    using RegionIdAndVariantId = std::pair<std::string, std::string>;
    const std::vector<RegionIdAndVariantId> getSortedIdPairs();

    const SampleParameters& sampleParams_;
    Reference& reference_;
    const RegionCatalog& regionCatalog_;
    const SampleFindings& sampleFindings_;
};

std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter);

}
