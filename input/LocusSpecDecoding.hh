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

#include <string>
#include <vector>

#include "boost/optional.hpp"

#include "common/GenomicRegion.hh"
#include "common/Parameters.hh"
#include "common/Reference.hh"
#include "region_spec/CNVLocusSpecification.hh"
#include "region_spec/GraphLocusSpecification.hh"

namespace ehunter
{

enum class VariantTypeFromUser
{
    kRareRepeat,
    kCommonRepeat,
    kSmallVariant,
    kSMN,
    kCNV
};

enum class VariantSubtypeFromUser
{
    kTarget,
    kBaseline
};

enum class LocusTypeFromUser
{
    kGraph,
    kCNV,
    kParalog
};

struct VariantDescriptionFromUser
{
    VariantDescriptionFromUser(
        std::string variantId, GenomicRegion variantLocation, VariantTypeFromUser variantType,
        boost::optional<VariantSubtypeFromUser> variantSubtype, boost::optional<std::string> variantStructure,
        boost::optional<bool> expectedNormalCN, boost::optional<double> regionGC,
        boost::optional<int> mappingQualityThreshold, boost::optional<int> maxCopyNumber,
        boost::optional<double> depthScaleFactor, boost::optional<double> standardDevidationOfCN2,
        boost::optional<std::vector<double>> meanDepthValues,
        boost::optional<std::vector<double>> priorCopyNumberFrequency)
        : variantId(std::move(variantId))
        , variantLocation(std::move(variantLocation))
        , variantType(std::move(variantType))
        , variantSubtype(variantSubtype)
        , variantStructure(variantStructure)
        , expectedNormalCN(expectedNormalCN)
        , regionGC(regionGC)
        , mappingQualityThreshold(mappingQualityThreshold)
        , maxCopyNumber(maxCopyNumber)
        , depthScaleFactor(depthScaleFactor)
        , standardDevidationOfCN2(standardDevidationOfCN2)
        , meanDepthValues(meanDepthValues)
        , priorCopyNumberFrequency(priorCopyNumberFrequency)

    {
    }

    std::string variantId;
    GenomicRegion variantLocation;
    VariantTypeFromUser variantType;
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
};

struct LocusDescriptionFromUser
{
    LocusDescriptionFromUser(
        std::string locusId, LocusTypeFromUser locusType, GenomicRegion locusLocation,
        std::vector<VariantDescriptionFromUser> variantDescriptionFromUsers, std::vector<GenomicRegion> targetRegions,
        std::vector<GenomicRegion> offtargetRegions, boost::optional<std::string> locusStructure,
        boost::optional<double> errorRate, boost::optional<double> likelihoodRatioThreshold,
        boost::optional<double> minLocusCoverage)
        : locusId(std::move(locusId))
        , locusType(std::move(locusType))
        , locusLocation(std::move(locusLocation))
        , variantDescriptionFromUsers(std::move(variantDescriptionFromUsers))
        , targetRegions(std::move(targetRegions))
        , offtargetRegions(std::move(offtargetRegions))
        , locusStructure(locusStructure)
        , errorRate(errorRate)
        , likelihoodRatioThreshold(likelihoodRatioThreshold)
        , minLocusCoverage(minLocusCoverage)

    {
    }

    std::string locusId;
    LocusTypeFromUser locusType;
    GenomicRegion locusLocation;
    std::vector<VariantDescriptionFromUser> variantDescriptionFromUsers;
    std::vector<GenomicRegion> targetRegions;
    std::vector<GenomicRegion> offtargetRegions;
    boost::optional<std::string> locusStructure;
    boost::optional<double> errorRate;
    boost::optional<double> likelihoodRatioThreshold;
    boost::optional<double> minLocusCoverage;
};

void assertValidity(const LocusDescriptionFromUser& userDescription);

GraphLocusSpecification
decodeGraphLocusSpecification(const LocusDescriptionFromUser& userDescription, const Reference& reference);

CNVLocusSpecification
decodeCNVLocusSpecification(const LocusDescriptionFromUser& userDescription, const Reference& reference);
}
