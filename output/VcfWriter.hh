//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>

#include "common/Parameters.hh"
#include "common/Reference.hh"
#include "locus_spec/CnvLocusSpec.hh"
#include "locus_spec/GraphLocusSpec.hh"
#include "locus_spec/LocusSpec.hh"
#include "workflow/LocusFindings.hh"

namespace ehunter
{

class GraphVariantVcfWriter : public VariantFindingsVisitor
{
public:
    GraphVariantVcfWriter(
        Reference& reference, const GraphLocusSpec& locusSpec, double locusDepth, const GraphVariantSpec& variantSpec,
        std::ostream& out)
        : reference_(reference)
        , locusSpec_(locusSpec)
        , locusDepth_(locusDepth)
        , variantSpec_(variantSpec)
        , out_(out)
    {
    }

    ~GraphVariantVcfWriter() = default;
    void visit(StrFindings& strFindings) override;
    void visit(SmallVariantFindings& smallVariantFindingsPtr) override;
    void visit(CnvVariantFindings& cnvVariantFindingsPtr) override;

private:
    Reference& reference_;
    const GraphLocusSpec& locusSpec_;
    double locusDepth_;
    const GraphVariantSpec& variantSpec_;
    std::ostream& out_;
};

/*
class CnvVariantVcfWriter : public VariantFindingsVisitor
{
public:
    CnvVariantVcfWriter(
        Reference& reference, const CnvLocusSpec& locusSpec, double locusDepth, std::ostream& out)
        : reference_(reference)
        , locusSpec_(locusSpec)
//, locusDepth_(locusDepth)
//, out_(out)
{
}

~CnvVariantVcfWriter() = default;
void visit(StrFindings& strFindings) override;
void visit(SmallVariantFindings& smallVariantFindingsPtr) override;
void visit(CnvVariantFindings& cnvVariantFindingsPtr) override;

private:
Reference& reference_;
const CnvLocusSpec& locusSpec_;
// double locusDepth_;
// std::ostream& out_;
};
*/

// TODO: Document the code after multi-unit repeat format is finalized (GT-598)
class VcfWriter
{
public:
    VcfWriter(
        std::string sampleId, Reference& reference, const LocusCatalog& regionCatalog,
        const SampleFindings& sampleFindings);

    friend std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter);

private:
    void writeHeader(std::ostream& out);
    void writeBody(std::ostream& out);
    using LocusIdAndVariantId = std::pair<std::string, std::string>;
    std::vector<LocusIdAndVariantId> getSortedIdPairs();

    std::string sampleId_;
    Reference& reference_;
    const LocusCatalog& regionCatalog_;
    const SampleFindings& sampleFindings_;
};

std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter);
}
