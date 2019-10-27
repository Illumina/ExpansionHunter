//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
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

#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"
#include "thirdparty/json/json.hpp"

#include "common/Common.hh"
#include "common/GenomicRegion.hh"
#include "common/Reference.hh"
#include "locus_spec/VariantSpecification.hh"

using std::make_shared;
using std::shared_ptr;

namespace ehunter
{

class LocusSpecification
{
public:
    LocusSpecification(
        std::string locusId, CopyNumberBySex copyNumberBySex, GenomicRegion locusLocation,
        GenotyperParameters genotyperParams)
        : locusId_(std::move(locusId))
        , copyNumberBySex_(copyNumberBySex)
        , locusLocation_(std::move(locusLocation))
        , parameters_(std::move(genotyperParams))
    {
    }

    virtual ~LocusSpecification() = default;

    const std::string& locusId() const { return locusId_; }
    CopyNumberBySex copyNumberBySex() const { return copyNumberBySex_; }
    const GenomicRegion& locusLocation() const { return locusLocation_; }
    const GenotyperParameters& genotyperParameters() const { return parameters_; }
    const std::vector<VariantSpecification>& variantSpecs() const { return variantSpecs_; }

    const VariantSpecification& getVariantSpecById(const std::string& variantSpecId) const;

protected:
    std::string locusId_;
    CopyNumberBySex copyNumberBySex_;
    GenomicRegion locusLocation_;
    GenotyperParameters parameters_;
    std::vector<VariantSpecification> variantSpecs_;
};
using LocusCatalog = std::map<std::string, shared_ptr<LocusSpecification>>;
}
