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

#include <boost/optional.hpp>

#include "common/GenomicRegion.hh"
#include "locus_spec/GraphLocusSpec.hh"

namespace ehunter
{

struct GraphVariantDecoding
{
    GraphVariantDecoding(std::string id, std::string type, GenomicRegion location)
        : id(std::move(id))
        , type(std::move(type))
        , location(location)
    {
    }
    std::string id;
    std::string type;
    GenomicRegion location;
};

struct GraphLocusDecoding
{
    std::string id;
    std::string structure;
    // std::string variantType;
    int flankLength;
    // Regions in the reference where we expect relevant reads to align
    std::vector<GenomicRegion> targetRegions;
    // Regions where additional relevant reads might be found that require filtering or special considerations
    std::vector<GenomicRegion> offtargetRegions;
    boost::optional<double> errorRate;
    boost::optional<double> likelihoodRatioThreshold;
    boost::optional<double> minLocusCoverage;

    std::vector<GraphVariantDecoding> variants;
};

GraphLocusSpec decode(const Reference& reference, const GraphLocusDecoding& locusEncoding);

}
