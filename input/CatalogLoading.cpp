//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "common/Common.hh"
#include "common/Reference.hh"
#include "common/WorkflowContext.hh"

using boost::optional;
using graphtools::NodeId;
using std::make_shared;
using std::map;
using std::ostream;
using std::shared_ptr;
using std::string;
using std::to_string;
using std::vector;

using Json = nlohmann::json;

namespace spd = spdlog;

namespace ehunter
{

enum class VariantSubtypeFromUser
{
    kTarget,
    kBaseline
};

enum class LocusTypeFromUser
{
    kGraph,
    kCNV,
    kParalog,
    kUnspecified
};

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

static Json makeArray(const Json& record)
{
    if (record.type() != Json::value_t::array)
    {
        return Json::array({ record });
    }
    else
    {
        return record;
    }
}

static LocusTypeFromUser decodeLocusTypeFromUser(const string& encoding)
{
    if (encoding == "Graph")
    {
        return LocusTypeFromUser::kGraph;
    }
    if (encoding == "CNV")
    {
        return LocusTypeFromUser::kCNV;
    }
    if (encoding == "Paralog")
    {
        return LocusTypeFromUser::kParalog;
    }
    else
    {
        throw std::logic_error("Encountered invalid locus type: " + encoding);
    }
}

static LocusTypeFromUser getLocusType(const Json& record)
{
    if (checkIfFieldExists(record, "LocusType"))
    {
        const auto encoding = record["LocusType"].get<string>();
        return decodeLocusTypeFromUser(encoding);
    }
    else
    {
        return LocusTypeFromUser::kUnspecified;
    }
}

static vector<string> generateIds(const string& locusId, const Json& variantLocationEncodings)
{
    if (variantLocationEncodings.size() == 1)
    {
        return { locusId };
    }

    vector<string> variantIds;
    for (const auto& locationEncoding : variantLocationEncodings)
    {
        string variantId = locusId;
        variantId += "_";
        variantId += locationEncoding.get<string>();
        variantIds.push_back(variantId);
    }

    return variantIds;
}

static GraphLocusEncoding loadGraphLocusEncoding(const Json& json, const Reference& reference)
{
    GraphLocusEncoding locus;

    assertFieldExists(json, "LocusId");
    locus.id = json["LocusId"].get<string>();

    assertFieldExists(json, "LocusStructure");
    locus.structure = json["LocusStructure"].get<string>();

    if (checkIfFieldExists(json, "TargetRegion"))
    {
        for (const auto& encoding : makeArray(json["TargetRegion"]))
        {
            auto region = decode(reference.contigInfo(), encoding.get<string>());
            locus.targetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(json, "OfftargetRegions"))
    {
        assertRecordIsArray(json["OfftargetRegions"]);
        for (const auto& encoding : json["OfftargetRegions"])
        {
            GenomicRegion region = decode(reference.contigInfo(), encoding.get<string>());
            locus.offtargetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(json, "ErrorRate"))
    {
        locus.errorRate = json["ErrorRate"].get<double>();
    }

    if (checkIfFieldExists(json, "LikelihoodRatioThreshold"))
    {
        locus.likelihoodRatioThreshold = json["LikelihoodRatioThreshold"].get<double>();
    }

    if (checkIfFieldExists(json, "MinimalLocusCoverage"))
    {
        locus.minLocusCoverage = json["MinimalLocusCoverage"].get<double>();
    }

    assertFieldExists(json, "Variants");
    makeArray(json["Variants"]);
    for (const auto& variant : json["Variants"])
    {
        assertFieldExists(variant, "ReferenceRegion");
        GenomicRegion region = decode(reference.contigInfo(), variant["ReferenceRegion"].get<string>());

        assertFieldExists(variant, "VariantType");
        auto variantType = variant["VariantType"].get<string>();

        std::string variantId;
        if (checkIfFieldExists(variant, "VariantId"))
        {
            variantId = variant["VariantId"].get<string>();
        }
        else
        {
            variantId = locus.id;
            variantId += "_";
            variantId += variant["ReferenceRegion"].get<string>();
        }

        locus.variants.emplace_back(GraphVariantEncoding(variantId, variantType, region));
    }
    return locus;
}

static CnvLocusEncoding loadCnvLocusEncoding(const Json& locusJson, const Reference& reference)
{
    CnvLocusEncoding cnvLocusDecoding;

    assertFieldExists(locusJson, "LocusId");
    auto locusId = locusJson["LocusId"].get<string>();
    cnvLocusDecoding.id = locusId;

    std::vector<CnvVariantEncoding> cnvAnalysisVariants;
    std::vector<CnvOutputVariantEncoding> cnvOutputVariants;

    assertFieldExists(locusJson, "OutputVariants");
    makeArray(locusJson["OutputVariants"]);
    for (const auto& variant : locusJson["OutputVariants"])
    {
        assertFieldExists(variant, "VariantId");
        auto variantId = variant["VariantId"].get<string>();

        assertFieldExists(variant, "ReferenceRegion");
        GenomicRegion region = decode(reference.contigInfo(), variant["ReferenceRegion"].get<string>());

        CnvOutputVariantEncoding cnvOutputVariantDecoding;
        cnvOutputVariantDecoding.id = variantId;
        cnvOutputVariantDecoding.location = region;
        cnvOutputVariants.emplace_back(cnvOutputVariantDecoding);
    }

    assertFieldExists(locusJson, "AnalysisVariants");
    makeArray(locusJson["AnalysisVariants"]);
    for (const auto& variant : locusJson["AnalysisVariants"])
    {
        assertFieldExists(variant, "ReferenceRegion");
        makeArray(variant["ReferenceRegion"]);
        std::vector<GenomicRegion> variantRegions;
        for (auto referenceRegion : variant["ReferenceRegion"])
        {
            GenomicRegion region = decode(reference.contigInfo(), referenceRegion.get<string>());
            variantRegions.push_back(region);
        }

        assertFieldExists(variant, "VariantId");
        string variantId = variant["VariantId"].get<string>();

        assertFieldExists(variant, "VariantSubtype");
        auto variantType = variant["VariantSubtype"].get<string>();

        bool expectedNormalCN;
        if (variantType == "Baseline")
        {
            assertFieldExists(variant, "ExpectedNormal");
            expectedNormalCN = variant["ExpectedNormal"].get<bool>();
        }
        else
        {
            expectedNormalCN = false;
        }

        assertFieldExists(variant, "GC");
        auto regionGC = variant["GC"].get<double>();

        assertFieldExists(variant, "MappingQualityThreshold");
        auto mappingQualityThreshold = variant["MappingQualityThreshold"].get<int>();

        assertFieldExists(variant, "MaxCopyNumber");
        auto maxCopyNumber = variant["MaxCopyNumber"].get<int>();

        assertFieldExists(variant, "DepthScaleFactor");
        auto depthScaleFactor = variant["DepthScaleFactor"].get<double>();

        assertFieldExists(variant, "StandardDeviationOfCN2");
        auto standardDeviationOfCN2 = variant["StandardDeviationOfCN2"].get<double>();

        assertFieldExists(variant, "MeanDepthValues");
        std::vector<double> meanDepths;
        for (const auto& encoding : variant["MeanDepthValues"])
        {
            meanDepths.push_back(encoding.get<double>());
        }
        auto meanDepthValues = meanDepths;

        assertFieldExists(variant, "PriorCopyNumberFreq");
        std::vector<double> priors;
        for (const auto& encoding : variant["PriorCopyNumberFreq"])
        {
            priors.push_back(encoding.get<double>());
        }
        auto priorCopyNumberFrequency = priors;

        CnvVariantEncoding variantDecoding;
        variantDecoding.id = variantId;
        variantDecoding.locations = variantRegions;
        variantDecoding.variantType = variantType;
        variantDecoding.expectedNormalCN = expectedNormalCN;
        variantDecoding.regionGC = regionGC;
        variantDecoding.mappingQualityThreshold = mappingQualityThreshold;
        variantDecoding.maxCopyNumber = maxCopyNumber;
        variantDecoding.depthScaleFactor = depthScaleFactor;
        variantDecoding.standardDeviationOfCN2 = standardDeviationOfCN2;
        variantDecoding.meanDepthValues = meanDepths;
        variantDecoding.priorCopyNumberFrequency = priorCopyNumberFrequency;
        cnvAnalysisVariants.emplace_back(variantDecoding);
    }
    cnvLocusDecoding.outputVariants = cnvOutputVariants;
    cnvLocusDecoding.variants = cnvAnalysisVariants;
    return cnvLocusDecoding;
}

static GraphLocusEncoding loadLegacyGraphLocusDecoding(const Json& record, const Reference& reference)
{
    GraphLocusEncoding locus;

    assertFieldExists(record, "LocusId");
    locus.id = record["LocusId"].get<string>();

    assertFieldExists(record, "LocusStructure");
    locus.structure = record["LocusStructure"].get<string>();

    if (checkIfFieldExists(record, "TargetRegion"))
    {
        for (const auto& encoding : makeArray(record["TargetRegion"]))
        {
            auto region = decode(reference.contigInfo(), encoding.get<string>());
            locus.targetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(record, "OfftargetRegions"))
    {
        assertRecordIsArray(record["OfftargetRegions"]);
        for (const auto& encoding : record["OfftargetRegions"])
        {
            GenomicRegion region = decode(reference.contigInfo(), encoding.get<string>());
            locus.offtargetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(record, "ErrorRate"))
    {
        locus.errorRate = record["ErrorRate"].get<double>();
    }

    if (checkIfFieldExists(record, "LikelihoodRatioThreshold"))
    {
        locus.likelihoodRatioThreshold = record["LikelihoodRatioThreshold"].get<double>();
    }

    if (checkIfFieldExists(record, "MinimalLocusCoverage"))
    {
        locus.minLocusCoverage = record["MinimalLocusCoverage"].get<double>();
    }

    vector<GenomicRegion> variantLocations;
    assertFieldExists(record, "ReferenceRegion");
    auto variantLocationJson = makeArray(record["ReferenceRegion"]);
    for (const auto& encoding : variantLocationJson)
    {
        GenomicRegion region = decode(reference.contigInfo(), encoding.get<string>());
        variantLocations.push_back(region);
    }

    vector<string> variantTypes;
    assertFieldExists(record, "VariantType");
    for (const auto& encoding : makeArray(record["VariantType"]))
    {
        variantTypes.push_back(encoding.get<string>());
    }

    if (variantTypes.size() != variantLocations.size())
    {
        throw std::runtime_error("Types and locations must be provided for each variant in locus " + locus.id);
    }

    vector<string> variantIds;
    if (checkIfFieldExists(record, "VariantId"))
    {
        for (const auto& variantId : makeArray(record["VariantId"]))
        {
            variantIds.push_back(variantId.get<string>());
        }
    }
    else
    {
        variantIds = generateIds(locus.id, variantLocationJson);
    }

    if (variantIds.size() != variantTypes.size())
    {
        throw std::runtime_error("An id must be provided for each variant in locus " + locus.id);
    }

    for (int index = 0; index != static_cast<int>(variantTypes.size()); ++index)
    {
        GraphVariantEncoding variant(variantIds[index], variantTypes[index], variantLocations[index]);
        locus.variants.push_back(variant);
    }

    return locus;
}

std::unique_ptr<LocusSpec> loadLocusSpec(const Json& record, const Reference& reference)
{
    LocusTypeFromUser locusType = getLocusType(record);
    switch (locusType)
    {
    case LocusTypeFromUser::kGraph:
        return decode(reference, loadGraphLocusEncoding(record, reference));
    case LocusTypeFromUser::kCNV:
        return decode(reference, loadCnvLocusEncoding(record, reference));
    case LocusTypeFromUser::kParalog:
        return decode(reference, loadCnvLocusEncoding(record, reference));
    default: // If locus type is not specified, assume that its legacy graph
        return decode(reference, loadLegacyGraphLocusDecoding(record, reference));
    }
}

LocusCatalog loadLocusCatalogFromDisk(const string& catalogPath, const Reference& reference)
{
    std::ifstream inputStream(catalogPath.c_str());

    if (!inputStream.is_open())
    {
        throw std::runtime_error("Failed to open catalog file " + catalogPath);
    }

    Json catalogJson;
    inputStream >> catalogJson;
    makeArray(catalogJson);

    WorkflowContext context;

    LocusCatalog catalog;
    for (const auto& locusJson : catalogJson)
    {
        try
        {
            std::unique_ptr<LocusSpec> locusSpec = loadLocusSpec(locusJson, reference);
            catalog.emplace(std::make_pair(locusSpec->locusId(), std::move(locusSpec)));
        }
        catch (const std::exception& except)
        {
            const string message = "Unable to load " + locusJson.dump() + ": " + except.what();
            if (context.heuristics().permissive())
            {
                spdlog::warn(message);
            }
            else
            {
                throw std::runtime_error(message);
            }
        }
    }

    return catalog;
}

std::vector<RegionInfo> loadNormRegionsFromDisk(const std::string& normRegionPath, const Reference& reference)
{
    std::vector<RegionInfo> normRegionInfo;
    std::ifstream inputStream(normRegionPath.c_str());

    if (!inputStream.is_open())
    {
        throw std::runtime_error("Failed to open norm region file " + normRegionPath);
    }

    Json normJson;
    inputStream >> normJson;
    makeArray(normJson);

    for (auto& regionJson : normJson)
    {
        assertFieldExists(regionJson, "GC");
        auto regionGC = regionJson["GC"].get<float>();
        assertFieldExists(regionJson, "ReferenceRegion");
        GenomicRegion region = decode(reference.contigInfo(), regionJson["ReferenceRegion"].get<string>());
        normRegionInfo.emplace_back(regionGC, region);
    }
    return normRegionInfo;
}
}
