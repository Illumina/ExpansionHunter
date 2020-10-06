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

#include "locus_spec/CnvLocusDecoding.hh"
#include "locus_spec/CnvLocusSpec.hh"

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

static CnvLocusType getCnvLocusType(const CnvLocusEncoding& encoding)
{
    CnvLocusType cnvLocusType = CnvLocusType::kNonoverlapping;
    for (const auto& variant : encoding.variants)
    {
        if (variant.variantType == "kBaseline" && !(variant.expectedNormalCN))
        {
            cnvLocusType = CnvLocusType::kOverlapping;
        }
    }
    return cnvLocusType;
}

static GenomicRegion getLocusLocation(const CnvLocusEncoding& locusEncoding)
{
    vector<GenomicRegion> variantLocations;
    for (const auto& variant : locusEncoding.variants)
    {
        for (auto region : *variant.locations)
        {
            variantLocations.push_back(region);
        }
    }
    return *variantLocations.begin();
}

static CnvVariantType getCnvVariantType(const CnvVariantEncoding variant)
{
    if (variant.variantType == "Baseline")
    {
        return CnvVariantType::kBaseline;
    }
    else if (variant.variantType == "Target")
    {
        return CnvVariantType::kTarget;
    }
    else
    {
        throw std::logic_error("Encountered invalid variant type: " + variant.variantType);
    }
}

std::unique_ptr<CnvLocusSpec> decode(const Reference& reference, const CnvLocusEncoding& encoding)
{
    GenomicRegion locusLocation = getLocusLocation(encoding);
    CopyNumberBySex copyNumberBySex = getCopyNumber(reference.contigInfo().getContigName(locusLocation.contigIndex()));
    CnvLocusType cnvLocusType = getCnvLocusType(encoding);

    CnvOutputVariant outputVariant;
    for (const auto& variant : encoding.outputVariants)
    {
        outputVariant.id = variant.id;
        outputVariant.location = variant.location;
    }

    unique_ptr<CnvLocusSpec> locusSpec(new CnvLocusSpec(encoding.id, cnvLocusType, copyNumberBySex, outputVariant));
    for (const auto& variant : encoding.variants)
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

        CnvVariantType variantType = getCnvVariantType(variant);
        locusSpec->addVariant(variant.id, variantType, *variant.locations, variantParameters);
    }

    return locusSpec;
}
}
