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

#include "common/Common.hh"
#include "common/GenomicRegion.hh"
#include "common/Parameters.hh"
#include "common/Reference.hh"
#include "locus_spec/LocusSpec.hh"
#include "locus_spec/CnvLocusSpec.hh"

namespace ehunter
{

enum class Base
{
    kA,
    kC,
    kG,
    kT
};

struct ParalogOutputVariant
{
    std::string id;
    boost::optional<GenomicRegion> location;
};

class SmallVariantSpec
{
public:
    SmallVariantSpec(
        std::string id, std::vector<GenomicRegion> locations, int mappingQualityThreshold, Base variantBase, Base nonvariantBase)
        : id_(std::move(id))
        , locations_(std::move(locations))
        , mappingQualityThreshold_(std::move(mappingQualityThreshold))
        , variantBase_(variantBase)
        , nonvariantBase_(nonvariantBase)
    {
        assertConsistency();
    }

    const std::string& id() const { return id_; }
    const std::vector<GenomicRegion>& locations() const { return locations_; }
    const int& mappingQualityThreshold() const { return mappingQualityThreshold_; }
    const Base& variantBase() const { return variantBase_; }
    const Base& nonvariantBase() const { return nonvariantBase_; }

    bool operator==(const SmallVariantSpec& other) const
    {
        return id_ == other.id_ && locations_ == other.locations_ && variantBase_ == other.variantBase_ && nonvariantBase_ == other.nonvariantBase_;
    }

    void assertConsistency() const;

private:
    std::string id_;
    std::vector<GenomicRegion> locations_;
    int mappingQualityThreshold_;
    Base variantBase_;
    Base nonvariantBase_;
};

class ParalogLocusSpec : public LocusSpec
{
public:
    ParalogLocusSpec(
        std::string locusId, CopyNumberBySex contigCopyNumber, ParalogOutputVariant outputVariant)
        : LocusSpec(locusId, contigCopyNumber)
        , outputVariant_(outputVariant)
    {
    }

    ~ParalogLocusSpec() override = default;

    std::vector<GenomicRegion> regionsWithReads() const override;
    const std::vector<CnvVariantSpec>& CnvVariants() const { return cnvVariants_; }
    const std::vector<SmallVariantSpec>& SmallVariants() const { return smallVariants_; }
    const ParalogOutputVariant& outputVariant() const { return outputVariant_; }
    void
    addCnvVariant(std::string id, CnvVariantType type, std::vector<GenomicRegion> referenceLocus, CnvGenotyperParameters parameters);
    void addSmallVariant(std::string id, std::vector<GenomicRegion> referenceLocus, int mappingQualityThreshold, Base variantBase, Base nonvariantBase);
    const GenomicRegion& getVariantLocationById(const std::string& id) const override;

private:
    std::vector<CnvVariantSpec> cnvVariants_;
    std::vector<SmallVariantSpec> smallVariants_;
    ParalogOutputVariant outputVariant_;
};
}
