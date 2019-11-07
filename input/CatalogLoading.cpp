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

static VariantSubtypeFromUser decodeVariantSubtypeFromUser(const string& encoding)
{
    if (encoding == "Target")
    {
        return VariantSubtypeFromUser::kTarget;
    }
    if (encoding == "Baseline")
    {
        return VariantSubtypeFromUser::kBaseline;
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

static CnvLocusDecoding loadCnvUserDescription(Json& locusJson, const ReferenceContigInfo& contigInfo)
{
    CnvLocusDecoding cnvLocusDecoding;

    assertFieldExists(locusJson, "LocusId");
    auto locusId = locusJson["LocusId"].get<string>();
    cnvLocusDecoding.id = locusId;

    std::vector<CnvVariantDecoding> cnvAnalysisVariants;
    std::vector<CnvOutputVariantDecoding> cnvOutputVariants;

    assertFieldExists(locusJson, "OutputVariants");
    makeArray(locusJson["OutputVariants"]);
    for (const auto& variant : locusJson["OutputVariants"])
    {
        assertFieldExists(variant, "VariantId");
        auto variantId = variant["VariantId"].get<string>();

        assertFieldExists(variant, "ReferenceRegion");
        GenomicRegion region = decode(contigInfo, variant["ReferenceRegion"].get<string>());

        CnvOutputVariantDecoding cnvOutputVariantDecoding = CnvOutputVariantDecoding(variantId, region);
        cnvOutputVariants.emplace_back(cnvOutputVariantDecoding);
    }

    assertFieldExists(locusJson, "AnalysisVariants");
    makeArray(locusJson["AnalysisVariants"]);
    for (const auto& variant : locusJson["AnalysisVariants"])
    {
        assertFieldExists(variant, "ReferenceRegion");
        GenomicRegion region = decode(contigInfo, variant["ReferenceRegion"].get<string>());

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
        auto standardDevidationOfCN2 = variant["StandardDeviationOfCN2"].get<double>();

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

        CnvVariantDecoding variantDecoding = CnvVariantDecoding(
            variantId, region, variantType, expectedNormalCN, regionGC, mappingQualityThreshold, maxCopyNumber,
            depthScaleFactor, standardDevidationOfCN2, meanDepths, priorCopyNumberFrequency);
        cnvAnalysisVariants.emplace_back(variantDecoding);
    }
    cnvLocusDecoding.outputVariants = cnvOutputVariants;
    cnvLocusDecoding.variants = cnvAnalysisVariants;
    return cnvLocusDecoding;
}

static LocusDescriptionFromUser loadUserDescription(Json& locusJson, const ReferenceContigInfo& contigInfo)
{
    assertFieldExists(locusJson, "LocusId");
    auto locusId = locusJson["LocusId"].get<string>();

    assertFieldExists(locusJson, "LocusType");
    LocusTypeFromUser locusType = decodeLocusTypeFromUser(locusJson["LocusType"].get<string>());

    boost::optional<std::string> locusStructure;
    if (locusType == LocusTypeFromUser::kGraph)
    {
        assertFieldExists(locusJson, "LocusStructure");
        locusStructure = locusJson["LocusStructure"].get<string>();
    }

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
            variantId = variant["VariantId"].get<string>();
        }
        else
        {
            variantId = locusId;
            variantId += "_";
            variantId += variant["ReferenceRegion"].get<string>();
        }

        boost::optional<VariantSubtypeFromUser> variantSubtype;
        boost::optional<std::string> variantStructure;
        boost::optional<bool> expectedNormalCN;
        boost::optional<double> regionGC;
        boost::optional<int> mappingQualityThreshold;
        boost::optional<int> maxCopyNumber;
        boost::optional<double> depthScaleFactor;
        boost::optional<double> standardDevidationOfCN2;
        boost::optional<std::vector<double>> meanDepthValues;
        boost::optional<std::vector<double>> priorCopyNumberFrequency;

        if (checkIfFieldExists(variant, "VariantStructure"))
        {
            variantStructure = variant["VariantStructure"].get<string>();
        }

        if (variantType == VariantTypeFromUser::kCNV)
        {
            assertFieldExists(variant, "VariantSubtype");
            variantSubtype = decodeVariantSubtypeFromUser(variant["VariantSubtype"].get<string>());

            if (variantSubtype == VariantSubtypeFromUser::kBaseline)
            {
                assertFieldExists(variant, "ExpectedNormal");
                expectedNormalCN = variant["ExpectedNormal"].get<bool>();
            }
            else
            {
                expectedNormalCN = false;
            }

            assertFieldExists(variant, "GC");
            regionGC = variant["GC"].get<double>();

            assertFieldExists(variant, "MappingQualityThreshold");
            mappingQualityThreshold = variant["MappingQualityThreshold"].get<int>();

            assertFieldExists(variant, "MaxCopyNumber");
            maxCopyNumber = variant["MaxCopyNumber"].get<int>();

            assertFieldExists(variant, "DepthScaleFactor");
            depthScaleFactor = variant["DepthScaleFactor"].get<double>();

            assertFieldExists(variant, "StandardDeviationOfCN2");
            standardDevidationOfCN2 = variant["StandardDeviationOfCN2"].get<double>();

            assertFieldExists(variant, "MeanDepthValues");
            std::vector<double> meanDepths;
            for (const auto& encoding : variant["MeanDepthValues"])
            {
                meanDepths.push_back(encoding.get<double>());
            }
            meanDepthValues = meanDepths;

            assertFieldExists(variant, "PriorCopyNumberFreq");
            std::vector<double> priors;
            for (const auto& encoding : variant["PriorCopyNumberFreq"])
            {
                priors.push_back(encoding.get<double>());
            }
            priorCopyNumberFrequency = priors;
        }
        VariantDescriptionFromUser variantDescription = VariantDescriptionFromUser(
            variantId, region, variantType, variantSubtype, variantStructure, expectedNormalCN, regionGC,
            mappingQualityThreshold, maxCopyNumber, depthScaleFactor, standardDevidationOfCN2, meanDepthValues,
            priorCopyNumberFrequency);
        variantDescriptions.push_back(variantDescription);
        variantLocations.push_back(region);
    }
    GenomicRegion locusLocation = mergeRegions(variantLocations, 50000);
    if (locusType == LocusTypeFromUser::kGraph)
    {
        locusLocation = mergeRegions(variantLocations, 500);
    }
    else if (locusType == LocusTypeFromUser::kCNV)
    {
        assertFieldExists(locusJson, "ReferenceRegion");
        locusLocation = decode(contigInfo, locusJson["ReferenceRegion"].get<string>());
    }

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
        locusId, locusType, locusLocation, variantDescriptions, targetRegions, offtargetRegions, locusStructure,
        errorRate, likelihoodRatioThreshold, minLocusCoverage);
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
    for (auto& locusJson : catalogJson)
    {
        // assertFieldExists(locusJson, "LocusType");
        // auto locusType = locusJson["LocusType"].get<string>();
        // if (locusType == "CNV")
        //{
        //    auto cnvUserDescription = loadCnvUserDescription(locusJson, reference.contigInfo());
        //    CnvLocusSpec cnvLocusSpec = decodeCnvLocus(reference, cnvUserDescription);
        //}
        LocusDescriptionFromUser userDescription = loadUserDescription(locusJson, reference.contigInfo());
        try
        {
            if (userDescription.locusType == LocusTypeFromUser::kGraph)
            {
                GraphLocusSpec locusSpec = decodeGraphLocusSpecification(userDescription, reference);
                catalog.emplace(std::make_pair(locusSpec.locusId(), make_shared<GraphLocusSpec>(locusSpec)));
            }
            // if (userDescription.locusType == LocusTypeFromUser::kCNV)
            //{
            //    CnvLocusSpec locusSpec = decodeCnvLocusSpecification(userDescription, reference);
            //    catalog.emplace(std::make_pair(locusSpec.locusId(), make_shared<CnvLocusSpec>(locusSpec)));
            //}
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
        normRegionInfo.push_back(RegionInfo(regionGC, region));
    }
    return normRegionInfo;
}
}
