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

#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "common/GenomicRegion.hh"
#include "locus_spec/CnvLocusSpec.hh"

namespace ehunter
{

struct CnvVariantDecoding
{
    std::string id;
    boost::optional<GenomicRegion> location;
    std::string variantType;
    bool expectedNormalCN;
    double regionGC;
    int mappingQualityThreshold;
    int maxCopyNumber;
    double depthScaleFactor;
    double standardDevidationOfCN2;
    std::vector<double> meanDepthValues;
    std::vector<double> priorCopyNumberFrequency;
};

struct CnvOutputVariantDecoding
{
    std::string id;
    boost::optional<GenomicRegion> location;
};

struct CnvLocusDecoding
{
    std::string id;
    std::vector<CnvVariantDecoding> variants;
    std::vector<CnvOutputVariantDecoding> outputVariants;
};

CnvLocusSpec decode(const Reference& reference, const CnvLocusDecoding& encoding);
}
