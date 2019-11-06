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

#include <stdexcept>

using std::runtime_error;
using std::string;
using std::to_string;
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

static CnvLocusType getCnvLocusType(const CnvLocusDecoding& encoding)
{
    CnvLocusType cnvLocusType = CnvLocusType::kNonoverlapping;
    for (const auto& variant : encoding.variants)
    {
        if (variant.variantType == "kBaseline" && !(*variant.expectedNormalCN))
        {
            cnvLocusType = CnvLocusType::kOverlapping;
        }
    }
    return cnvLocusType;
}

CnvLocusSpec decode(const Reference& reference, const CnvLocusDecoding& encoding)
{
    auto copyNumberBySex = getCopyNumber(reference.contigInfo().getContigName(locusLocation.contigIndex()));
    auto cnvLocusType = getCnvLocusType(encoding);

    CnvLocusSpec locusSpec(encoding.id, cnvLocusType, copyNumberBySex);
    for (const auto& variant : encoding.variants)
    {
        auto parameters = CnvGenotyperParameters(
            variant.regionGC, variant.maxCopyNumber, variant.mappingQualityThreshold, variant.depthScaleFactor,
            variant.standardDevidationOfCN2, variant.meanDepthValues, variant.priorCopyNumberFrequency,
            variant.expectedNormalCN);
        locusSpec.addVariant(variant.id, variant.location, parameters);
    }

    return locusSpec;
}
}

