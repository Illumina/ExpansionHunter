//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "input/CatalogLoading.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "thirdparty/json/json.hpp"
#include "thirdparty/spdlog/spdlog.h"

#include "common/Common.hh"
#include "common/Reference.hh"
#include "input/GraphBlueprint.hh"
#include "input/RegionGraph.hh"

using boost::optional;
using graphtools::NodeId;
using std::map;
using std::ostream;
using std::string;
using std::to_string;
using std::vector;

using Json = nlohmann::json;

namespace spd = spdlog;

namespace ehunter
{

enum class VariantDescriptionFromUser
{
    kRareRepeat,
    kCommonRepeat,
    kSmallVariant,
    kSMN
};

enum class InputRecordType
{
    kRegionWithSingleRepeat,
    kRegionWithMultipleRepeats,
    kUnknown
};

AlleleCount determineExpectedAlleleCount(Sex sex, const string& chrom)
{
    const bool isFemaleChromY = sex == Sex::kFemale && (chrom == "chrY" || chrom == "Y");
    if (isFemaleChromY)
    {
        return AlleleCount::kZero;
    }

    const bool isSexChrom = chrom == "chrX" || chrom == "X" || chrom == "chrY" || chrom == "Y";
    if (sex == Sex::kMale && isSexChrom)
    {
        return AlleleCount::kOne;
    }

    return AlleleCount::kTwo;
}

static bool checkIfFieldExists(const Json& record, const string& fieldName)
{
    return record.find(fieldName) != record.end();
}

static void assertFieldExists(const Json& record, const string& fieldName)
{
    if (!checkIfFieldExists(record, fieldName))
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Field " + fieldName + " must be present in " + out.str());
    }
}

static void assertRecordIsArray(const Json& record)
{
    if (!record.is_array())
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Expected array but got this instead " + out.str());
    }
}

static VariantDescriptionFromUser decodeVariantDescription(const string& encoding)
{
    if (encoding == "RareRepeat")
    {
        return VariantDescriptionFromUser::kRareRepeat;
    }
    if (encoding == "Repeat")
    {
        return VariantDescriptionFromUser::kCommonRepeat;
    }
    if (encoding == "SmallVariant")
    {
        return VariantDescriptionFromUser::kSmallVariant;
    }
    if (encoding == "SMN")
    {
        return VariantDescriptionFromUser::kSMN;
    }
    else
    {
        throw std::logic_error("Encountered invalid variant type: " + encoding);
    }
}

static vector<string> combine(const std::string& prefix, const vector<string>& suffixes)
{
    vector<string> combinedStrings;
    for (const auto& suffix : suffixes)
    {
        combinedStrings.push_back(prefix + "_" + suffix);
    }

    return combinedStrings;
}

static void makeArray(Json& record)
{
    if (record.type() != Json::value_t::array)
    {
        record = Json::array({ record });
    }
}

static bool doesFeatureDefineVariant(GraphBlueprintFeatureType featureType)
{
    switch (featureType)
    {
    case GraphBlueprintFeatureType::kInsertionOrDeletion:
    case GraphBlueprintFeatureType::kSkippableRepeat:
    case GraphBlueprintFeatureType::kUnskippableRepeat:
    case GraphBlueprintFeatureType::kSwap:
        return true;

    case GraphBlueprintFeatureType::kLeftFlank:
    case GraphBlueprintFeatureType::kRightFlank:
    case GraphBlueprintFeatureType::kInterruption:
        return false;

    default:
        std::stringstream encoding;
        encoding << featureType;
        throw std::logic_error("Unrecognized feature type: " + encoding.str());
    }
}

static std::size_t countVariants(const GraphBlueprint& blueprint)
{
    std::size_t numVariants = 0;
    for (const auto& feature : blueprint)
    {
        if (doesFeatureDefineVariant(feature.type))
        {
            ++numVariants;
        }
    }

    return numVariants;
}

static VariantType determineVariantType(GraphBlueprintFeatureType featureType)
{
    switch (featureType)
    {
    case GraphBlueprintFeatureType::kInsertionOrDeletion:
    case GraphBlueprintFeatureType::kSwap:
        return VariantType::kSmallVariant;
    case GraphBlueprintFeatureType::kSkippableRepeat:
    case GraphBlueprintFeatureType::kUnskippableRepeat:
        return VariantType::kRepeat;
    default:
        std::ostringstream encoding;
        encoding << featureType;
        throw std::logic_error("Feature of type " + encoding.str() + " does not define a variant");
    }
}

static VariantSubtype determineVariantSubtype(
    GraphBlueprintFeatureType featureType, VariantDescriptionFromUser userDescription, const Region referenceRegion)
{
    if (featureType == GraphBlueprintFeatureType::kInsertionOrDeletion)
    {
        if (referenceRegion.length() == 0)
        {
            return VariantSubtype::kInsertion;
        }
        else
        {
            return VariantSubtype::kDeletion;
        }
    }
    else if (featureType == GraphBlueprintFeatureType::kSwap)
    {
        if (userDescription == VariantDescriptionFromUser::kSMN)
        {
            return VariantSubtype::kSMN;
        }
        else
        {
            return VariantSubtype::kSwap;
        }
    }
    else if (userDescription == VariantDescriptionFromUser::kCommonRepeat)
    {
        return VariantSubtype::kCommonRepeat;
    }
    else if (userDescription == VariantDescriptionFromUser::kRareRepeat)
    {
        return VariantSubtype::kRareRepeat;
    }
    else
    {
        std::ostringstream message;
        message << featureType;
        throw std::logic_error("Feature " + message.str() + " does not correspond to variant");
    }
}


static optional<NodeId>
determineReferenceNode(const GraphBlueprintFeature& feature, const Reference& reference, const Region& referenceRegion)
{
    const string refSequence
        = reference.getSequence(referenceRegion.chrom(), referenceRegion.start(), referenceRegion.end());

    optional<NodeId> optionalReferenceNode;
    for (int index = 0; index != static_cast<int>(feature.nodeIds.size()); ++index)
    {
        if (refSequence == feature.sequences[index])
        {
            optionalReferenceNode = feature.nodeIds[index];
            break;
        }
    }

    return optionalReferenceNode;
}

static GraphBlueprint generateBlueprint(const Reference& reference, const Region& region, const string& locusStructure)
{
    // Reference repeat flanks should be at least as long as reads.
    const int kFlankLen = 1500;
    const int64_t leftFlankStart = region.start() - kFlankLen;
    const int64_t rightFlankEnd = region.end() + kFlankLen;

    const string leftFlank = reference.getSequence(region.chrom(), leftFlankStart, region.start());
    const string regionSequence = reference.getSequence(region.chrom(), region.start(), region.end());
    const string rightFlank = reference.getSequence(region.chrom(), region.end(), rightFlankEnd);

    return decodeFeaturesFromRegex(leftFlank + locusStructure + rightFlank);
}

static Region mergeRegions(const vector<Region>& regions)
{
    const int kMaxMergeDistance = 500;
    vector<Region> mergedReferenceRegions = merge(regions, kMaxMergeDistance);
    if (mergedReferenceRegions.size() != 1)
    {
        std::stringstream out;
        for (const Region& region : regions)
        {
            out << region << " ";
        }
        throw std::runtime_error(
            "Expected reference regions to be closer than " + std::to_string(kMaxMergeDistance)
            + " from one another: " + out.str());
    }

    return mergedReferenceRegions.front();
}

static LocusSpecification loadLocusSpecification(Json& locusJson, Sex sampleSex, const Reference& reference)
{
    assertFieldExists(locusJson, "LocusId");
    const string& locusId = locusJson["LocusId"];

    assertFieldExists(locusJson, "ReferenceRegion");
    makeArray(locusJson["ReferenceRegion"]);
    vector<Region> referenceRegions;
    for (const string& encoding : locusJson["ReferenceRegion"])
    {
        referenceRegions.emplace_back(encoding);
    }

    vector<string> variantIds = combine(locusId, locusJson["ReferenceRegion"]);

    assertFieldExists(locusJson, "LocusStructure");
    const string& locusStructure = locusJson["LocusStructure"];

    assertFieldExists(locusJson, "VariantType");
    makeArray(locusJson["VariantType"]);
    vector<VariantDescriptionFromUser> variantDescriptions;
    for (const string& encoding : locusJson["VariantType"])
    {
        variantDescriptions.push_back(decodeVariantDescription(encoding));
    }

    const Region mergedReferenceRegion = mergeRegions(referenceRegions);
    vector<Region> targetRegions;
    if (checkIfFieldExists(locusJson, "TargetRegion"))
    {
        makeArray(locusJson["TargetRegion"]);
        for (const string& locusEncoding : locusJson["TargetRegion"])
        {
            targetRegions.emplace_back(locusEncoding);
        }
    }
    else
    {
        targetRegions = { mergedReferenceRegion };
    }

    vector<Region> offtargetRegions;
    if (checkIfFieldExists(locusJson, "OfftargetRegions"))
    {
        assertRecordIsArray(locusJson["OfftargetRegions"]);
        for (const string& locusEncoding : locusJson["OfftargetRegions"])
        {
            offtargetRegions.push_back(Region(locusEncoding));
        }
    }

    GraphBlueprint blueprint = generateBlueprint(reference, mergedReferenceRegion, locusStructure);
    graphtools::Graph locusGraph = makeRegionGraph(blueprint);

    std::size_t numVariants = countVariants(blueprint);
    if (numVariants == 0)
    {
        std::stringstream out;
        out << locusJson;
        throw std::runtime_error("Locus must contain at least one variant: " + out.str());
    }

    if (referenceRegions.size() != numVariants || variantDescriptions.size() != numVariants)
    {
        std::stringstream out;
        out << locusJson;
        throw std::runtime_error("Expected reference region and type for each variant: " + out.str());
    }

    AlleleCount expectedAlleleCount = determineExpectedAlleleCount(sampleSex, mergedReferenceRegion.chrom());
    LocusSpecification regionSpec(locusId, targetRegions, expectedAlleleCount, locusGraph);

    int variantIndex = 0;
    for (const auto& feature : blueprint)
    {
        if (doesFeatureDefineVariant(feature.type))
        {
            const Region& referenceRegion = referenceRegions[variantIndex];
            VariantDescriptionFromUser variantDescription = variantDescriptions[variantIndex];

            VariantType variantType = determineVariantType(feature.type);
            VariantSubtype variantSubtype = determineVariantSubtype(feature.type, variantDescription, referenceRegion);
            optional<NodeId> optionalReferenceNode = determineReferenceNode(feature, reference, referenceRegion);

            VariantClassification classification(variantType, variantSubtype);
            regionSpec.addVariantSpecification(
                variantIds[variantIndex], classification, referenceRegion, feature.nodeIds, optionalReferenceNode);
            ++variantIndex;
        }
    }

    return regionSpec;
}

RegionCatalog loadRegionCatalogFromDisk(const string& catalogPath, const Reference& reference, Sex sampleSex)
{
    std::ifstream inputStream(catalogPath.c_str());

    if (!inputStream.is_open())
    {
        throw std::runtime_error("Failed to open catalog file " + catalogPath);
    }

    Json catalogJson;
    inputStream >> catalogJson;
    makeArray(catalogJson);

    RegionCatalog catalog;
    for (auto& locusJson : catalogJson)
    {
        LocusSpecification locusSpec = loadLocusSpecification(locusJson, sampleSex, reference);
        catalog.emplace(std::make_pair(locusSpec.regionId(), locusSpec));
    }

    return catalog;
}

}
