//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Xiao Chen <xchen2@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "locus_spec/ParalogLocusSpec.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "common/Common.hh"
#include "common/Reference.hh"

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

vector<GenomicRegion> ParalogLocusSpec::regionsWithReads() const
{
    vector<GenomicRegion> regions;
    for (const auto& variant : cnvVariants_)
    {
        for (const auto& region : variant.locations())
        {
            regions.push_back(region);
        }
    }
    /*
    for (const auto& variant : smallVariants_)
    {
        regions.push_back(variant.locations().geneALocation);
        regions.push_back(variant.locations().geneBLocation);
    }
    */

    return regions;
}

void ParalogLocusSpec::addCnvVariant(
    std::string id, CnvVariantType type, std::vector<GenomicRegion> referenceLocus, CnvGenotyperParameters parameters)
{
    cnvVariants_.emplace_back(std::move(id), type, std::move(referenceLocus), std::move(parameters));
}

void ParalogLocusSpec::addSmallVariant(
    std::string id, std::vector<GenomicRegion> referenceLocus, int mappingQualityThreshold, std::pair<Base, Base> bases)
{
    smallVariants_.emplace_back(
        std::move(id), SmallVariantLocations(*referenceLocus.begin(), *referenceLocus.end()), mappingQualityThreshold,
        SmallVariantBases(bases.first, bases.second));
}

void SmallVariantSpec::assertConsistency() const
{
    bool variantIsValid = (bases_.geneABase != bases_.geneBBase);

    if (!variantIsValid)
    {
        throw std::logic_error("Definition of variant " + id_ + " is inconsistent");
    }
}

const GenomicRegion& ParalogLocusSpec::getVariantLocationById(const string& id) const
{
    if (outputVariant_.id == id)
    {
        return *(outputVariant_.location);
    }

    throw std::logic_error("There is no variant " + id + " in locus " + locusId_);
}
}
