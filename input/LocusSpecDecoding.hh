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
#include "region_spec/LocusSpecification.hh"

namespace ehunter
{

enum class VariantTypeFromUser
{
    kRareRepeat,
    kCommonRepeat,
    kSmallVariant,
    kSMN
};

struct LocusDescriptionFromUser
{
    LocusDescriptionFromUser(
        std::string locusId, std::string locusStructure, std::vector<std::string> variantIds,
        GenomicRegion locusLocation, std::vector<GenomicRegion> variantLocations,
        std::vector<GenomicRegion> targetRegions, std::vector<GenomicRegion> offtargetRegions,
        std::vector<VariantTypeFromUser> variantTypesFromUser, boost::optional<double> errorRate,
        boost::optional<double> likelihoodRatioThreshold, boost::optional<double> minLocusCoverage)
        : locusId(std::move(locusId))
        , locusStructure(std::move(locusStructure))
        , variantIds(std::move(variantIds))
        , locusLocation(std::move(locusLocation))
        , variantLocations(std::move(variantLocations))
        , targetRegions(std::move(targetRegions))
        , offtargetRegions(std::move(offtargetRegions))
        , variantTypesFromUser(std::move(variantTypesFromUser))
        , errorRate(errorRate)
        , likelihoodRatioThreshold(likelihoodRatioThreshold)
        , minLocusCoverage(minLocusCoverage)

    {
    }

    std::string locusId;
    std::string locusStructure;
    std::vector<std::string> variantIds;
    GenomicRegion locusLocation;
    std::vector<GenomicRegion> variantLocations;
    std::vector<GenomicRegion> targetRegions;
    std::vector<GenomicRegion> offtargetRegions;
    std::vector<VariantTypeFromUser> variantTypesFromUser;
    boost::optional<double> errorRate;
    boost::optional<double> likelihoodRatioThreshold;
    boost::optional<double> minLocusCoverage;
};

void assertValidity(const LocusDescriptionFromUser& userDescription);

LocusSpecification
decodeLocusSpecification(const LocusDescriptionFromUser& userDescription, const Reference& reference);

}
