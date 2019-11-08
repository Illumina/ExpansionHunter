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

namespace ehunter
{

enum class CnvLocusType
{
    kOverlapping,
    kNonoverlapping
};

enum class CnvVariantType
{
    kTarget,
    kBaseline
};

// Per-variant parameters for CNV variant genotyping
struct CnvGenotyperParameters
{
    double regionGC;
    int maxCopyNumber;
    int mappingQualityThreshold;
    double depthScaleFactor;
    double standardDeviationOfCN2;
    std::vector<double> meanDepthValues;
    std::vector<double> priorCopyNumberFrequency;
    bool expectedNormal;
};

struct CnvOutputVariant
{
    std::string id;
    boost::optional<GenomicRegion> location;
};

class CnvVariantSpec
{
public:
    CnvVariantSpec(std::string id, CnvVariantType variantType, GenomicRegion location, CnvGenotyperParameters genotyperParams)
        : id_(std::move(id))
        , variantType_(std::move(variantType))
        , location_(std::move(location))
        , genotyperParams_(std::move(genotyperParams))
    {
        assertConsistency();
    }

    const std::string& id() const { return id_; }
    const GenomicRegion& location() const { return location_; }
    const CnvVariantType& variantType() const { return variantType_; }
    const CnvGenotyperParameters& genotyperParams() const { return genotyperParams_; }

    bool operator==(const CnvVariantSpec& other) const 
    { 
        return id_ == other.id_ && variantType_ == other.variantType_ && location_ == other.location_; 
    }

    void assertConsistency() const;

private:
    std::string id_;
    CnvVariantType variantType_;
    GenomicRegion location_;
    CnvGenotyperParameters genotyperParams_;
};

class CnvLocusSpec : public LocusSpec
{
public:
    CnvLocusSpec(
        std::string locusId, CnvLocusType locusType, CopyNumberBySex contigCopyNumber, CnvOutputVariant outputVariant)
        : LocusSpec(locusId, contigCopyNumber)
        , locusType_(locusType)
        , outputVariant_(outputVariant)
    {
    }

    ~CnvLocusSpec() override = default;

    std::vector<GenomicRegion> regionsWithReads() const;
    const CnvLocusType& locusType() const { return locusType_; }
    const std::vector<CnvVariantSpec>& variants() const { return variants_; }
    const CnvOutputVariant& outputVariant() const { return outputVariant_; }
    void addVariant(std::string id, CnvVariantType type, GenomicRegion referenceLocus, CnvGenotyperParameters parameters);

private:
    CnvLocusType locusType_;
    std::vector<CnvVariantSpec> variants_;
    CnvOutputVariant outputVariant_;
};
}
