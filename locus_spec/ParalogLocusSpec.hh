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
#include "locus_spec/CnvLocusSpec.hh"
#include "locus_spec/LocusSpec.hh"

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

struct SmallVariantLocations
{
    SmallVariantLocations(GenomicRegion geneALocation, GenomicRegion geneBLocation)
        : geneALocation(geneALocation)
        , geneBLocation(geneBLocation)
    {
    }

    bool operator==(const SmallVariantLocations& other) const
    {
        return (other.geneALocation == geneALocation && other.geneBLocation == geneBLocation);
    }

    GenomicRegion geneALocation;
    GenomicRegion geneBLocation;
};

struct SmallVariantBases
{
    SmallVariantBases(Base geneABase, Base geneBBase)
        : geneABase(geneABase)
        , geneBBase(geneBBase)
    {
    }

    bool operator==(const SmallVariantBases& other) const
    {
        return (other.geneABase == geneABase && other.geneBBase == geneBBase);
    }

    Base geneABase;
    Base geneBBase;
};

class SmallVariantSpec
{
public:
    SmallVariantSpec(
        std::string id, SmallVariantLocations locations, int mappingQualityThreshold, SmallVariantBases bases)
        : id_(std::move(id))
        , locations_(std::move(locations))
        , mappingQualityThreshold_(std::move(mappingQualityThreshold))
        , bases_(bases)
    {
        assertConsistency();
    }

    const std::string& id() const { return id_; }
    const SmallVariantLocations& locations() const { return locations_; }
    const int& mappingQualityThreshold() const { return mappingQualityThreshold_; }
    const SmallVariantBases& variantBases() const { return bases_; }

    bool operator==(const SmallVariantSpec& other) const
    {
        return id_ == other.id_ && locations_ == other.locations_ && bases_ == other.bases_;
    }

    void assertConsistency() const;

private:
    std::string id_;
    SmallVariantLocations locations_;
    int mappingQualityThreshold_;
    SmallVariantBases bases_;
};

class ParalogLocusSpec : public LocusSpec
{
public:
    ParalogLocusSpec(std::string locusId, CopyNumberBySex contigCopyNumber, ParalogOutputVariant outputVariant)
        : LocusSpec(locusId, contigCopyNumber)
        , outputVariant_(outputVariant)
    {
    }

    ~ParalogLocusSpec() override = default;

    std::vector<GenomicRegion> regionsWithReads() const override;
    const std::vector<CnvVariantSpec>& cnvVariants() const { return cnvVariants_; }
    const std::vector<SmallVariantSpec>& smallVariants() const { return smallVariants_; }
    const ParalogOutputVariant& outputVariant() const { return outputVariant_; }
    void addCnvVariant(
        std::string id, CnvVariantType type, std::vector<GenomicRegion> referenceLocus,
        CnvGenotyperParameters parameters);
    void addSmallVariant(
        std::string id, std::vector<GenomicRegion> referenceLocus, int mappingQualityThreshold,
        std::pair<Base, Base> bases);
    const GenomicRegion& getVariantLocationById(const std::string& id) const override;

private:
    std::vector<CnvVariantSpec> cnvVariants_;
    std::vector<SmallVariantSpec> smallVariants_;
    ParalogOutputVariant outputVariant_;
};
}
