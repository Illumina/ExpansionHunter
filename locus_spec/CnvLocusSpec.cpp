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

#include "locus_spec/CnvLocusSpec.hh"

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

vector<GenomicRegion> CnvLocusSpec::regionsWithReads() const
{
    vector<GenomicRegion> regions;
    for (const auto& variant : variants_)
    {
        regions.push_back(variant.location());
    }

    return regions;
}

void CnvLocusSpec::addVariant(
    std::string id, CnvVariantType type, GenomicRegion referenceLocus, CnvGenotyperParameters parameters)
{
    variants_.emplace_back(std::move(id), type, std::move(referenceLocus), std::move(parameters));
}

void CnvVariantSpec::assertConsistency() const
{
    bool variantIsValid = (variantType_ == CnvVariantType::kBaseline || variantType_ == CnvVariantType::kTarget);

    if (!variantIsValid)
    {
        throw std::logic_error("Definition of variant " + id_ + " is inconsistent");
    }
}

const GenomicRegion& CnvLocusSpec::getVariantLocationById(const string& id) const
{
    if (outputVariant_.id == id)
    {
        return *(outputVariant_.location);
    }

    throw std::logic_error("There is no variant " + id + " in locus " + locusId_);
}
}
