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

#include "core/GenomicRegion.hh"
#include "core/Parameters.hh"
#include "core/Reference.hh"
#include "locus/LocusSpecification.hh"

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
    std::string locusId;

    /// Regular expression defining the structure of the locus
    std::string locusStructure;
    std::vector<std::string> variantIds;
    std::vector<GenomicRegion> referenceRegions;
    std::vector<GenomicRegion> targetRegions;
    std::vector<GenomicRegion> offtargetRegions;
    std::vector<VariantTypeFromUser> variantTypesFromUser;
    boost::optional<double> errorRate;
    boost::optional<double> likelihoodRatioThreshold;
    boost::optional<double> minLocusCoverage;

    /// If true, turn on additional motif processing for a repeat
    bool useRFC1MotifAnalysis = false;
};

void assertValidity(const LocusDescriptionFromUser& userDescription);

LocusSpecification decodeLocusSpecification(
    const LocusDescriptionFromUser& userDescription, const Reference& reference,
    const HeuristicParameters& heuristicParams);

}
