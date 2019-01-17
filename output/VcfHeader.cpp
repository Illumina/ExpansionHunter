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

#include "output/VcfHeader.hh"

#include <memory>

using std::ostream;
using std::string;

namespace ehunter
{

void FieldDescriptionWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    if (!repeatFindingsPtr->optionalGenotype())
    {
        return;
    }

    tryAddingFieldDescription(FieldType::kInfo, "SVTYPE", "1", "String", "Type of structural variant");
    tryAddingFieldDescription(FieldType::kInfo, "END", "1", "Integer", "End position of the variant");
    tryAddingFieldDescription(FieldType::kInfo, "REF", "1", "Integer", "Reference copy number");
    tryAddingFieldDescription(FieldType::kInfo, "RL", "1", "Integer", "Reference length in bp");
    tryAddingFieldDescription(FieldType::kInfo, "RU", "1", "String", "Repeat unit in the reference orientation");

    const string kRepidFieldDescription = "Repeat identifier from the repeat specification file";
    tryAddingFieldDescription(FieldType::kInfo, "REPID", "1", "String", kRepidFieldDescription);
    tryAddingFieldDescription(FieldType::kFormat, "GT", "1", "String", "Genotype");

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

    tryAddingFieldDescription(FieldType::kFilter, "PASS", "", "", "All filters passed");

    const auto repeatNodeId = variantSpec_.nodes().front();
    const string& repeatUnit = regionSpec_.regionGraph().nodeSeq(repeatNodeId);
    const auto& referenceLocus = variantSpec_.referenceLocus();
    const int referenceSize = referenceLocus.length() / repeatUnit.length();

    const RepeatGenotype& genotype = repeatFindingsPtr->optionalGenotype().get();

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

void FieldDescriptionWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    if (!smallVariantFindingsPtr->optionalGenotype())
    {
        return;
    }

    tryAddingFieldDescription(FieldType::kFormat, "GT", "1", "String", "Genotype");
    tryAddingFieldDescription(FieldType::kFilter, "PASS", "", "", "All filters passed");
}

void FieldDescriptionWriter::tryAddingFieldDescription(
    FieldType fieldType, const string& id, const string& number, const string& contentType, const string& description)
{
    const auto key = std::make_pair(fieldType, id);
    if (fieldDescriptions_.find(key) == fieldDescriptions_.end())
    {
        FieldDescription fieldDescription(
            fieldType, id, std::move(number), std::move(contentType), std::move(description));
        fieldDescriptions_.emplace(std::make_pair(std::move(key), std::move(fieldDescription)));
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

void outputVcfHeader(const RegionCatalog& regionCatalog, const SampleFindings& sampleFindings, ostream& out)
{
    out << "##fileformat=VCFv4.1\n";

    FieldDescriptionCatalog fieldDescriptionCatalog;

    for (const auto& regionIdAndFindings : sampleFindings)
    {
        const string& regionId = regionIdAndFindings.first;
        const LocusSpecification& regionSpec = regionCatalog.at(regionId);
        const RegionFindings& regionFindings = regionIdAndFindings.second;

        for (const auto& variantIdAndFindings : regionFindings)
        {
            const string& variantId = variantIdAndFindings.first;
            const VariantSpecification& variantSpec = regionSpec.getVariantSpecById(variantId);

            FieldDescriptionWriter descriptionWriter(regionSpec, variantSpec);
            variantIdAndFindings.second->accept(&descriptionWriter);
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
