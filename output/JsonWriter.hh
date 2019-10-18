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

#include "common/Parameters.hh"
#include "region_spec/GraphLocusSpecification.hh"
#include "region_spec/CNVLocusSpecification.hh"
#include "region_spec/LocusSpecification.hh"
#include "workflow/LocusFindings.hh"

#include "thirdparty/json/json.hpp"

namespace ehunter
{

class GraphVariantJsonWriter : public VariantFindingsVisitor
{
public:
    GraphVariantJsonWriter(
        const ReferenceContigInfo& contigInfo, const GraphLocusSpecification& locusSpec,
        const VariantSpecification& variantSpec)
        : contigInfo_(contigInfo)
        , locusSpec_(locusSpec)
        , variantSpec_(variantSpec)
    {
    }

    ~GraphVariantJsonWriter() = default;
    void visit(StrFindings& strFindings) override;
    void visit(SmallVariantFindings& smallVariantFindings) override;
    void visit(CNVVariantFindings& cnvVariantFindings) override;
    nlohmann::json record() const { return record_; }

private:
    const ReferenceContigInfo& contigInfo_;
    const GraphLocusSpecification& locusSpec_;
    const VariantSpecification& variantSpec_;
    nlohmann::json record_;
};


class CNVVariantJsonWriter : public VariantFindingsVisitor
{
public:
    CNVVariantJsonWriter(
        const ReferenceContigInfo& contigInfo, const CNVLocusSpecification& locusSpec)
        : contigInfo_(contigInfo)
        , locusSpec_(locusSpec)
    {
    }

    ~CNVVariantJsonWriter() = default;
    void visit(StrFindings& strFindings) override;
    void visit(SmallVariantFindings& smallVariantFindings) override;
    void visit(CNVVariantFindings& cnvVariantFindings) override;
    nlohmann::json record() const { return record_; }

private:
    const ReferenceContigInfo& contigInfo_;
    const CNVLocusSpecification& locusSpec_;
    nlohmann::json record_;
};


class JsonWriter
{
public:
    JsonWriter(
        const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
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
