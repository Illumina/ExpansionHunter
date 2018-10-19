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

#include "stats/ReadSupportCalculator.hh"

using std::map;
using std::string;
using Json = nlohmann::json;
using std::to_string;

std::ostream& operator<<(std::ostream& out, JsonWriter& jsonWriter)
{
    jsonWriter.write(out);
    return out;
}

JsonWriter::JsonWriter(
    const string& sampleName, int readLength, const RegionCatalog& regionSpecs, const SampleFindings& sampleFindings)
    : sampleName_(sampleName)
    , readLength_(readLength)
    , regionSpecs_(regionSpecs)
    , sampleFindings_(sampleFindings)
{
}

void JsonWriter::write(std::ostream& out)
{
    Json array = Json::array();
    for (const auto& regionIdAndRegionFindings : sampleFindings_)
    {
        const string& regionId = regionIdAndRegionFindings.first;
        const RegionFindings& regionFindings = regionIdAndRegionFindings.second;

        addRegionFindings(regionId, regionFindings, array);
    };

    out << std::setw(4) << array << std::endl;
}

void JsonWriter::addRegionFindings(const std::string& regionId, const RegionFindings& regionFindings, Json& array)
{
    const RegionSpec& regionSpec = regionSpecs_.at(regionId);

    for (const auto& blueprintComponent : regionSpec.regionBlueprint())
    {
        if (blueprintComponent.type() == RegionBlueprintComponent::Type::kRepeat)
        {
            const auto& repeatIdAndRepeatFindingsIter = regionFindings.find(blueprintComponent.id());
            assert(repeatIdAndRepeatFindingsIter != regionFindings.end());
            const RepeatFindings& repeatFindings = repeatIdAndRepeatFindingsIter->second;
            if (repeatFindings.optionalGenotype())
            {
                addRepeatFindings(blueprintComponent, repeatFindings, array);
            }
        }
    }
}

template <typename T> static string streamToString(const T& streamableObject)
{
    std::stringstream out;
    out << streamableObject;
    return out.str();
}

string encodeGenotype(const RepeatGenotype& genotype)
{
    string encoding = std::to_string(genotype.shortAlleleSizeInUnits());

    if (genotype.numAlleles() == 2)
    {
        encoding = encoding + "/" + std::to_string(genotype.longAlleleSizeInUnits());
    }

    return encoding;
}

static string encodeSupportCounts(const ReadSupportCalculator& readSupportCalculator, int repeatSize)
{
    const string countsOfSpanningReads = to_string(readSupportCalculator.getCountOfConsistentSpanningReads(repeatSize));
    const string countsOfFlankingReads = to_string(readSupportCalculator.getCountOfConsistentFlankingReads(repeatSize));
    return countsOfSpanningReads + "-" + countsOfFlankingReads;
}

string JsonWriter::encodeRepeatAlleleSupport(const string& repeatUnit, const RepeatFindings& repeatFindings)
{
    const int maxUnitsInRead = std::ceil(readLength_ / static_cast<double>(repeatUnit.length()));
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

void JsonWriter::addRepeatFindings(
    const RegionBlueprintComponent& repeatBlueprint, const RepeatFindings& repeatFindings, Json& array)
{
    Json record;

    const auto referenceRegion = repeatBlueprint.referenceRegion();
    assert(referenceRegion);

    record["RepeatId"] = repeatBlueprint.id();
    if (repeatBlueprint.referenceRegion())
    {
        record["ReferenceLocus"] = streamToString(*repeatBlueprint.referenceRegion());
    }

    const string repeatUnit = repeatBlueprint.sequence();
    record["RepeatUnit"] = repeatUnit;
    assert(repeatFindings.optionalGenotype());
    record["Genotype"] = encodeGenotype(*repeatFindings.optionalGenotype());
    record["GenotypeConfidenceInterval"] = streamToString(*repeatFindings.optionalGenotype());
    record["GenotypeSupport"] = encodeRepeatAlleleSupport(repeatUnit, repeatFindings);
    record["CountsOfSpanningReads"] = streamToString(repeatFindings.countsOfSpanningReads());
    record["CountsOfFlankingReads"] = streamToString(repeatFindings.countsOfFlankingReads());

    array.push_back(record);
}
