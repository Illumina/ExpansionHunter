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
#include "input/LocusSpecDecoding.hh"

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

static void makeArray(Json& record)
{
    if (record.type() != Json::value_t::array)
    {
        record = Json::array({ record });
    }
}

static VariantTypeFromUser decodeVariantTypeFromUser(const string& encoding)
{
    if (encoding == "RareRepeat")
    {
        return VariantTypeFromUser::kRareRepeat;
    }
    if (encoding == "Repeat")
    {
        return VariantTypeFromUser::kCommonRepeat;
    }
    if (encoding == "SmallVariant")
    {
        return VariantTypeFromUser::kSmallVariant;
    }
    if (encoding == "SMN")
    {
        return VariantTypeFromUser::kSMN;
    }
    if (encoding == "CNV")
    {
        return VariantTypeFromUser::kCNV;
    }
    else
    {
        throw std::logic_error("Encountered invalid variant type: " + encoding);
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

/*
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
*/
static GenomicRegion mergeRegions(const vector<GenomicRegion>& regions)
{
    const int kMaxMergeDistance = 500;
    vector<GenomicRegion> mergedReferenceRegions = merge(regions, kMaxMergeDistance);
    if (mergedReferenceRegions.size() != 1)
    {
        std::stringstream out;
        for (const GenomicRegion& region : regions)
        {
            out << region << " ";
        }
        throw std::runtime_error(
            "Expected reference regions to be closer than " + to_string(kMaxMergeDistance)
            + " from one another: " + out.str());
    }

    return mergedReferenceRegions.front();
}

static LocusDescriptionFromUser loadUserDescription(Json& locusJson, const ReferenceContigInfo& contigInfo)
{
    assertFieldExists(locusJson, "LocusId");
    auto locusId = locusJson["LocusId"].get<string>();

    assertFieldExists(locusJson, "LocusType");
    LocusTypeFromUser locusType = decodeLocusTypeFromUser(locusJson["LocusType"].get<string>());

    assertFieldExists(locusJson, "LocusStructure");
    auto locusStructure = locusJson["LocusStructure"].get<string>();

    vector<GenomicRegion> variantLocations;

    assertFieldExists(locusJson, "Variants");
    makeArray(locusJson["Variants"]);
    vector<VariantDescriptionFromUser> variantDescriptions;
    for (const auto& variant : locusJson["Variants"])
    {
        assertFieldExists(variant, "VariantType");
        auto variantType = decodeVariantTypeFromUser(variant["VariantType"].get<string>());

        assertFieldExists(variant, "ReferenceRegion");
        GenomicRegion region = decode(contigInfo, variant["ReferenceRegion"].get<string>());

        string variantId;
        if (checkIfFieldExists(variant, "VariantId"))
        {
            variantId = variant["VariantId"];
        }
        else
        {
            variantId = locusId;
            variantId += "_";
            variantId += variant["ReferenceRegion"].get<string>();
        }

        VariantDescriptionFromUser variantDescription = VariantDescriptionFromUser(variantId, region, variantType);
        variantDescriptions.push_back(variantDescription);
        variantLocations.push_back(region);
    }

    GenomicRegion locusLocation = mergeRegions(variantLocations);

    vector<GenomicRegion> targetRegions;
    if (checkIfFieldExists(locusJson, "TargetRegion"))
    {
        makeArray(locusJson["TargetRegion"]);
        for (const auto& encoding : locusJson["TargetRegion"])
        {
            GenomicRegion region = decode(contigInfo, encoding.get<string>());
            targetRegions.push_back(region);
        }
    }

    vector<GenomicRegion> offtargetRegions;
    if (checkIfFieldExists(locusJson, "OfftargetRegions"))
    {
        assertRecordIsArray(locusJson["OfftargetRegions"]);
        for (const auto& encoding : locusJson["OfftargetRegions"])
        {
            GenomicRegion region = decode(contigInfo, encoding.get<string>());
            offtargetRegions.push_back(region);
        }
    }

    boost::optional<double> errorRate;
    if (checkIfFieldExists(locusJson, "ErrorRate"))
    {
        errorRate = locusJson["ErrorRate"].get<double>();
    }

    boost::optional<double> likelihoodRatioThreshold;
    if (checkIfFieldExists(locusJson, "LikelihoodRatioThreshold"))
    {
        likelihoodRatioThreshold = locusJson["LikelihoodRatioThreshold"].get<double>();
    }

    boost::optional<double> minLocusCoverage;
    if (checkIfFieldExists(locusJson, "MinimalLocusCoverage"))
    {
        minLocusCoverage = locusJson["MinimalLocusCoverage"].get<double>();
    }

    return LocusDescriptionFromUser(
        locusId, locusType, locusStructure, locusLocation, variantDescriptions, targetRegions, offtargetRegions,
        errorRate, likelihoodRatioThreshold, minLocusCoverage);
}

RegionCatalog loadLocusCatalogFromDisk(const string& catalogPath, const Reference& reference)
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

    RegionCatalog catalog;
    for (auto& locusJson : catalogJson)
    {
        LocusDescriptionFromUser userDescription = loadUserDescription(locusJson, reference.contigInfo());
        try
        {
            LocusSpecification locusSpec = decodeLocusSpecification(userDescription, reference);
            catalog.emplace(std::make_pair(locusSpec.locusId(), locusSpec));
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
}
