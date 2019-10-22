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

#include "LocusSpecification.hh"
#include "common/Common.hh"
#include "common/GenomicRegion.hh"
#include "common/Reference.hh"
#include "region_spec/VariantSpecification.hh"

namespace ehunter
{

class CnvLocusSpecification : public LocusSpecification
{
public:
    CnvLocusSpecification(
        RegionId locusId, LocusType locusType, CnvLocusSubtype locusSubtype, ContigCopyNumber contigCopyNumber,
        GenomicRegion locusLocation, GenotyperParameters genotyperParams)
        : LocusSpecification(locusId, locusType, contigCopyNumber, locusLocation, genotyperParams)
        , locusSubtype_(std::move(locusSubtype))
    {
    }
    ~CnvLocusSpecification() override = default;

    const CnvLocusSubtype& locusSubtype() const { return locusSubtype_; }
    void addVariantSpecification(
        std::string id, VariantClassification classification, GenomicRegion referenceLocus,
        boost::optional<CnvGenotyperParameters> paramters);

private:
    CnvLocusSubtype locusSubtype_;
};
}
