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
#include <string>

#include "locus_spec/GraphLocusSpec.hh"
#include "locus_spec/LocusSpec.hh"
#include "workflow/LocusFindings.hh"

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
    FieldDescriptionWriter(std::shared_ptr<LocusSpec> locusSpecPtr)
        : locusSpecPtr_(std::move(locusSpecPtr))
    {
    }

    ~FieldDescriptionWriter() = default;

    void visit(const StrFindings& strFindings) override;
    void visit(const SmallVariantFindings& findings) override;
    void visit(const CnvVariantFindings& findings) override;
    void visit(const ParalogSmallVariantFindings& findings) override;

    void tryAddingFieldDescription(
        FieldType fieldType, const std::string& id, const std::string& number, const std::string& contentType,
        const std::string& description);

    void dumpTo(FieldDescriptionCatalog& descriptionCatalog);

private:
    void addCommonFields();

    std::shared_ptr<LocusSpec> locusSpecPtr_;
    FieldDescriptionCatalog fieldDescriptions_;
};

void outputVcfHeader(const LocusCatalog& regionCatalog, const SampleFindings& sampleFindings, std::ostream& out);

std::ostream& operator<<(std::ostream& out, FieldType fieldType);
std::ostream& operator<<(std::ostream& out, const FieldDescription& fieldDescription);
}
