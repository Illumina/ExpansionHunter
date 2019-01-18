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

#include "output/JsonWriter.hh"

#include <sstream>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/optional.hpp>

#include "stats/ReadSupportCalculator.hh"

namespace ehunter
{

using std::map;
using std::string;
using Json = nlohmann::json;
using boost::optional;
using std::to_string;
using std::vector;

std::ostream& operator<<(std::ostream& out, JsonWriter& jsonWriter)
{
    jsonWriter.write(out);
    return out;
}

JsonWriter::JsonWriter(
    const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
    const SampleFindings& sampleFindings)
    : sampleParams_(sampleParams)
    , contigInfo_(contigInfo)
    , regionCatalog_(regionCatalog)
    , sampleFindings_(sampleFindings)
{
}

void JsonWriter::write(std::ostream& out)
{
    Json array = Json::array();

    for (const auto& regionIdAndFindings : sampleFindings_)
    {
        const string& regionId = regionIdAndFindings.first;
        const LocusSpecification& regionSpec = regionCatalog_.at(regionId);
        const RegionFindings& regionFindings = regionIdAndFindings.second;

        for (const auto& variantIdAndFindings : regionFindings)
        {
            const string& variantId = variantIdAndFindings.first;
            const VariantSpecification& variantSpec = regionSpec.getVariantSpecById(variantId);

            VariantJsonWriter variantWriter(sampleParams_, contigInfo_, regionSpec, variantSpec);
            variantIdAndFindings.second->accept(&variantWriter);
            array.push_back(variantWriter.record());
        }
    };

    out << std::setw(4) << array << std::endl;
}

template <typename T> static string streamToString(const T& streamableObject)
{
    std::stringstream out;
    out << streamableObject;
    return out.str();
}

static string encodeSupportCounts(const ReadSupportCalculator& readSupportCalculator, int repeatSize)
{
    const string countsOfSpanningReads = to_string(readSupportCalculator.getCountOfConsistentSpanningReads(repeatSize));
    const string countsOfFlankingReads = to_string(readSupportCalculator.getCountOfConsistentFlankingReads(repeatSize));
    return countsOfSpanningReads + "-" + countsOfFlankingReads;
}

static string encodeRepeatAlleleSupport(const string& repeatUnit, const RepeatFindings& repeatFindings, int readLength)
{
    const int maxUnitsInRead = std::ceil(readLength / static_cast<double>(repeatUnit.length()));
    ReadSupportCalculator readSupportCalculator(
        maxUnitsInRead, repeatFindings.countsOfSpanningReads(), repeatFindings.countsOfFlankingReads());

    assert(repeatFindings.optionalGenotype());
    const RepeatGenotype& genotype = *repeatFindings.optionalGenotype();

    string supportEncoding = encodeSupportCounts(readSupportCalculator, genotype.shortAlleleSizeInUnits());

    if (genotype.numAlleles() == 2)
    {
        supportEncoding += "/" + encodeSupportCounts(readSupportCalculator, genotype.longAlleleSizeInUnits());
    }

    return supportEncoding;
}

static string encodeGenotype(const RepeatGenotype& genotype)
{
    string encoding = std::to_string(genotype.shortAlleleSizeInUnits());

    if (genotype.numAlleles() == 2)
    {
        encoding = encoding + "/" + std::to_string(genotype.longAlleleSizeInUnits());
    }

    return encoding;
}

void VariantJsonWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    assert(variantSpec_.classification().type == VariantType::kRepeat);

    const RepeatFindings& repeatFindings = *repeatFindingsPtr;

    record_.clear();
    record_["VariantId"] = variantSpec_.id();
    record_["ReferenceRegion"] = encode(contigInfo_, variantSpec_.referenceLocus());
    record_["VariantType"] = streamToString(variantSpec_.classification().type);
    record_["VariantSubtype"] = streamToString(variantSpec_.classification().subtype);

    const auto repeatNodeId = variantSpec_.nodes().front();
    const auto& repeatUnit = regionSpec_.regionGraph().nodeSeq(repeatNodeId);
    record_["RepeatUnit"] = repeatUnit;

    record_["CountsOfSpanningReads"] = streamToString(repeatFindings.countsOfSpanningReads());
    record_["CountsOfFlankingReads"] = streamToString(repeatFindings.countsOfFlankingReads());
    record_["CountsOfInrepeatReads"] = streamToString(repeatFindings.countsOfInrepeatReads());

    if (repeatFindings.optionalGenotype() != boost::none)
    {
        record_["Genotype"] = encodeGenotype(*repeatFindings.optionalGenotype());
        record_["GenotypeConfidenceInterval"] = streamToString(*repeatFindings.optionalGenotype());
        record_["GenotypeSupport"] = encodeRepeatAlleleSupport(repeatUnit, repeatFindings, sampleParams_.readLength());
    }
}

void VariantJsonWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    const SmallVariantFindings& indelFindings = *smallVariantFindingsPtr;
    record_.clear();
    record_["VariantId"] = variantSpec_.id();
    record_["VariantType"] = streamToString(variantSpec_.classification().type);
    record_["VariantSubtype"] = streamToString(variantSpec_.classification().subtype);
    record_["ReferenceRegion"] = streamToString(variantSpec_.referenceLocus());
    record_["CountOfRefReads"] = indelFindings.numRefReads();
    record_["CountOfAltReads"] = indelFindings.numAltReads();
    record_["StatusOfRefAllele"] = streamToString(indelFindings.refAllelePresenceStatus());
    record_["StatusOfAltAllele"] = streamToString(indelFindings.altAllelePresenceStatus());
    if (indelFindings.optionalGenotype() != boost::none)
    {
        record_["Genotype"] = streamToString(*indelFindings.optionalGenotype());
    }
}

}
