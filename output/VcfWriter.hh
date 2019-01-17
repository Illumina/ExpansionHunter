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

#include "common/Reference.hh"
#include "region_analysis/VariantFindings.hh"
#include "region_spec/LocusSpecification.hh"

namespace ehunter
{

class VariantVcfWriter : public VariantFindingsVisitor
{
public:
    VariantVcfWriter(
        const LocusSpecification& regionSpec, const VariantSpecification& variantSpec, int readLength,
        Reference& reference, std::ostream& out)
        : regionSpec_(regionSpec)
        , variantSpec_(variantSpec)
        , readLength_(readLength)
        , reference_(reference)
        , out_(out)
    {
    }

    ~VariantVcfWriter() = default;
    void visit(const RepeatFindings* repeatFindingsPtr) override;
    void visit(const SmallVariantFindings* smallVariantFindingsPtr) override;

private:
    const LocusSpecification& regionSpec_;
    const VariantSpecification& variantSpec_;
    const int readLength_;
    Reference& reference_;
    std::ostream& out_;
};

// TODO: Document the code after multi-unit repeat format is finalized (GT-598)
class VcfWriter
{
public:
    VcfWriter(
        const std::string& sampleName, int readLength, const RegionCatalog& regionCatalog,
        const SampleFindings& sampleFindings, Reference& reference);

    friend std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter);

private:
    void writeHeader(std::ostream& out);
    void writeBody(std::ostream& out);

    const std::string sampleName_;
    const int readLength_;
    const RegionCatalog& regionCatalog_;
    const SampleFindings& sampleFindings_;
    Reference& reference_;
};

std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter);

}
