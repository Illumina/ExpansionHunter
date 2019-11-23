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

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include "locus_spec/ParalogLocusDecoding.hh"
#include "locus_spec/ParalogLocusSpec.hh"

#include <stdexcept>

using std::runtime_error;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;

namespace ehunter
{

static CopyNumberBySex getCopyNumber(const string& contig)
{
    if (contig == "chrY" || contig == "Y")
    {
        return CopyNumberBySex::kZeroInFemaleOneInMale;
    }

    if (contig == "chrX" || contig == "X")
    {
        return CopyNumberBySex::kTwoInFemaleOneInMale;
    }

    return CopyNumberBySex::kTwoInFemaleTwoInMale;
}

static GenomicRegion getLocusLocation(const ParalogLocusEncoding& locusEncoding)
{
    vector<GenomicRegion> variantLocations;
    for (const auto& variant : locusEncoding.cnvVariants)
    {
        for (auto region : *variant.locations)
        {
            variantLocations.push_back(region);
        }
    }
    return *variantLocations.begin();
}

static Base decodeBase(const std::string base)
{
    if (base == "A")
    {
        return Base::kA;
    }
    else if (base == "C")
    {
        return Base::kC;
    }
    else if (base == "T")
    {
        return Base::kT;
    }
    else if (base == "G")
    {
        return Base::kG;
    }
    
    throw std::logic_error("Variant base " + base + " is not recognized.");
}

static std::pair<Base, Base> getSmallVariantBases(const std::string variantStructure)
{
    vector<string> components;
    boost::algorithm::split(components, variantStructure, boost::algorithm::is_any_of("()|"));

    if (components.size() != 2)
    {
        throw std::logic_error("Unexpected small variant structure format: " + variantStructure);
    }

    Base variantBase = decodeBase(components[0]);
    Base nonvariantBase = decodeBase(components[1]);
    return std::pair<Base, Base>(variantBase, nonvariantBase);
}

std::unique_ptr<ParalogLocusSpec> decode(const Reference& reference, const ParalogLocusEncoding& encoding)
{
    GenomicRegion locusLocation = getLocusLocation(encoding);
    CopyNumberBySex copyNumberBySex = getCopyNumber(reference.contigInfo().getContigName(locusLocation.contigIndex()));

    ParalogOutputVariant outputVariant;
    for (const auto& variant : encoding.outputVariants)
    {
        outputVariant.id = variant.id;
        outputVariant.location = variant.location;
    }

    unique_ptr<ParalogLocusSpec> locusSpec(new ParalogLocusSpec(encoding.id, copyNumberBySex, outputVariant));
    for (const auto& variant : encoding.cnvVariants)
    {
        CnvGenotyperParameters variantParameters;
        variantParameters.regionGC = variant.regionGC;
        variantParameters.mappingQualityThreshold = variant.mappingQualityThreshold;
        variantParameters.maxCopyNumber = variant.maxCopyNumber;
        variantParameters.depthScaleFactor = variant.depthScaleFactor;
        variantParameters.standardDeviationOfCN2 = variant.standardDeviationOfCN2;
        variantParameters.meanDepthValues = variant.meanDepthValues;
        variantParameters.priorCopyNumberFrequency = variant.priorCopyNumberFrequency;
        variantParameters.expectedNormal = variant.expectedNormalCN;

        CnvVariantType variantType = CnvVariantType::kTarget;
        locusSpec->addCnvVariant(variant.id, variantType, *variant.locations, variantParameters);
    }
    for (const auto& variant : encoding.smallVariants)
    {
        auto variantStructure = variant.variantStructure;
        auto variantBases = getSmallVariantBases(variantStructure);
        locusSpec->addSmallVariant(variant.id, *variant.locations, variant.mappingQualityThreshold, variantBases.first, variantBases.second);
    }

    return locusSpec;
}
}
