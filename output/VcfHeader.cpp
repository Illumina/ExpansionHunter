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

#include "output/VcfHeader.hh"

#include <memory>

using std::ostream;
using std::string;

namespace ehunter
{

void FieldDescriptionWriter::addCommonFields()
{
    const string kVaridFieldDescription = "Variant identifier as specified in the variant catalog";
    tryAddingFieldDescription(FieldType::kInfo, "VARID", "1", "String", kVaridFieldDescription);
    tryAddingFieldDescription(FieldType::kFormat, "GT", "1", "String", "Genotype");
    tryAddingFieldDescription(FieldType::kFormat, "LC", "1", "Float", "Locus coverage");
    tryAddingFieldDescription(FieldType::kFilter, "PASS", "", "", "All filters passed");
}

void FieldDescriptionWriter::visit(const StrFindings& strFindings)
{
    if (!strFindings.optionalGenotype())
    {
        return;
    }

    auto graphLocusSpecPtr = std::static_pointer_cast<GraphLocusSpec>(locusSpecPtr_);
    const GraphVariantSpec& variantSpec = graphLocusSpecPtr->getVariantById(strFindings.variantId());

    addCommonFields();
    tryAddingFieldDescription(FieldType::kInfo, "SVTYPE", "1", "String", "Type of structural variant");
    tryAddingFieldDescription(FieldType::kInfo, "END", "1", "Integer", "End position of the variant");
    tryAddingFieldDescription(FieldType::kInfo, "REF", "1", "Integer", "Reference copy number");
    tryAddingFieldDescription(FieldType::kInfo, "RL", "1", "Integer", "Reference length in bp");
    tryAddingFieldDescription(FieldType::kInfo, "RU", "1", "String", "Repeat unit in the reference orientation");

    const string kRepidFieldDescription = "Repeat identifier as specified in the variant catalog";
    tryAddingFieldDescription(FieldType::kInfo, "REPID", "1", "String", kRepidFieldDescription);

    const string kSoFieldDescription
        = "Type of reads that support the allele; can be SPANNING, FLANKING, or INREPEAT meaning that the reads span, "
          "flank, or are fully contained in the repeat";
    tryAddingFieldDescription(FieldType::kFormat, "SO", "1", "String", kSoFieldDescription);

    const string kRepcnFieldDescription = "Number of repeat units spanned by the allele";
    tryAddingFieldDescription(FieldType::kFormat, "REPCN", "1", "String", kRepcnFieldDescription);
    tryAddingFieldDescription(FieldType::kFormat, "REPCI", "1", "String", "Confidence interval for REPCN");

    const string kAdflFieldDescription = "Number of flanking reads consistent with the allele";
    tryAddingFieldDescription(FieldType::kFormat, "ADFL", "1", "String", kAdflFieldDescription);

    const string kAdspFieldDescription = "Number of spanning reads consistent with the allele";
    tryAddingFieldDescription(FieldType::kFormat, "ADSP", "1", "String", kAdspFieldDescription);

    const string kAdirFieldDescription = "Number of in-repeat reads consistent with the allele";
    tryAddingFieldDescription(FieldType::kFormat, "ADIR", "1", "String", kAdirFieldDescription);

    const auto repeatNodeId = variantSpec.nodes().front();
    const string& repeatUnit = graphLocusSpecPtr->graph().nodeSeq(repeatNodeId);
    const auto& referenceLocus = variantSpec.location();
    const int referenceSize = referenceLocus.length() / repeatUnit.length();

    const RepeatGenotype& genotype = strFindings.optionalGenotype().get();

    if (genotype.shortAlleleSizeInUnits() != referenceSize)
    {
        const string sizeEncoding = std::to_string(genotype.shortAlleleSizeInUnits());
        const string description = "Allele comprised of " + sizeEncoding + " repeat units";
        tryAddingFieldDescription(FieldType::kAlt, "STR" + sizeEncoding, "", "", description);
    }

    if (genotype.longAlleleSizeInUnits() != referenceSize)
    {
        const string sizeEncoding = std::to_string(genotype.longAlleleSizeInUnits());
        const string description = "Allele comprised of " + sizeEncoding + " repeat units";
        tryAddingFieldDescription(FieldType::kAlt, "STR" + sizeEncoding, "", "", description);
    }
}

void FieldDescriptionWriter::visit(const ParalogSmallVariantFindings& findings)
{
    auto paralogFindings = findings;
}

void FieldDescriptionWriter::visit(const CnvVariantFindings& cnvFindings)
{
    auto finding = cnvFindings;
    tryAddingFieldDescription(
        FieldType::kInfo, "VARID", "1", "String", "Variant identifier as specified in the variant catalog");
    tryAddingFieldDescription(FieldType::kFilter, "PASS", "", "", "All filters passed");
    tryAddingFieldDescription(FieldType::kFilter, "LowQ", "", "", "CNV call is Low quality");
    tryAddingFieldDescription(FieldType::kFilter, "BaselineFail", "", "", "CNV call failed check against baseline regions");
    tryAddingFieldDescription(FieldType::kInfo, "CNLEN", "1", "Integer", "Number of reference positions spanned by this CNV");
    tryAddingFieldDescription(FieldType::kFormat, "CN", "1", "Integer", "Absolute Copy number");
    tryAddingFieldDescription(FieldType::kFormat, "CNC", "1", "Integer", "Copy number change");
}

void FieldDescriptionWriter::visit(const SmallVariantFindings& findings)
{
    if (!findings.optionalGenotype())
    {
        return;
    }

    auto graphLocusSpecPtr = std::static_pointer_cast<GraphLocusSpec>(locusSpecPtr_);
    const auto& variantSpec = graphLocusSpecPtr->getVariantById(findings.variantId());

    addCommonFields();
    tryAddingFieldDescription(
        FieldType::kFormat, "AD", ".", "Integer", "Allelic depths for the ref and alt alleles in the order listed");
    if (variantSpec.classification().subtype == GraphVariantClassification::Subtype::kSMN)
    {
        tryAddingFieldDescription(
            FieldType::kFormat, "RPL", "1", "Float", "Log-Likelihood ratio for the presence of the reference allele");
        tryAddingFieldDescription(
            FieldType::kFormat, "DST", "1", "Character",
            "Result ('+' detected, '-' undetected, '?' undetermined) of the test represented by the variant");
    }
}

void FieldDescriptionWriter::tryAddingFieldDescription(
    FieldType fieldType, const string& id, const string& number, const string& contentType, const string& description)
{
    const auto key = std::make_pair(fieldType, id);
    if (fieldDescriptions_.find(key) == fieldDescriptions_.end())
    {
        FieldDescription fieldDescription(fieldType, id, number, contentType, description);
        fieldDescriptions_.emplace(std::make_pair(key, std::move(fieldDescription)));
    }
}

void FieldDescriptionWriter::dumpTo(FieldDescriptionCatalog& descriptionCatalog)
{
    descriptionCatalog.insert(fieldDescriptions_.begin(), fieldDescriptions_.end());
}

FieldDescription::FieldDescription(
    FieldType fieldType, string id, string number, string contentType, string description)
    : fieldType(fieldType)
    , id(std::move(id))
    , number(std::move(number))
    , contentType(std::move(contentType))
    , description(std::move(description))
{
}

void outputVcfHeader(const LocusCatalog& locusCatalog, const SampleFindings& sampleFindings, ostream& out)
{
    out << "##fileformat=VCFv4.1\n";

    FieldDescriptionCatalog fieldDescriptionCatalog;

    for (const auto& locusIdAndFindings : sampleFindings)
    {
        const string& locusId = locusIdAndFindings.first;
        auto LocusSpecPtr = locusCatalog.at(locusId);

        const LocusFindings& locusFindings = locusIdAndFindings.second;

        for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
        {
            FieldDescriptionWriter descriptionWriter(LocusSpecPtr);
            variantIdAndFindings.second->accept(descriptionWriter);
            descriptionWriter.dumpTo(fieldDescriptionCatalog);
        }
    }

    for (const auto& fieldIdAndDescription : fieldDescriptionCatalog)
    {
        const auto& description = fieldIdAndDescription.second;
        out << description << "\n";
    }
}

std::ostream& operator<<(std::ostream& out, FieldType fieldType)
{
    switch (fieldType)
    {
    case FieldType::kAlt:
        out << "ALT";
        break;
    case FieldType::kFormat:
        out << "FORMAT";
        break;
    case FieldType::kInfo:
        out << "INFO";
        break;
    case FieldType::kFilter:
        out << "FILTER";
        break;
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, const FieldDescription& fieldDescription)
{
    switch (fieldDescription.fieldType)
    {
    case FieldType::kInfo:
    case FieldType::kFormat:
        out << "##" << fieldDescription.fieldType << "=<ID=" << fieldDescription.id
            << ",Number=" << fieldDescription.number << ",Type=" << fieldDescription.contentType << ",Description=\""
            << fieldDescription.description << "\">";
        break;
    case FieldType::kAlt:
    case FieldType::kFilter:
        out << "##" << fieldDescription.fieldType << "=<ID=" << fieldDescription.id << ",Description=\""
            << fieldDescription.description << "\">";
        break;
    }

    return out;
}
}
