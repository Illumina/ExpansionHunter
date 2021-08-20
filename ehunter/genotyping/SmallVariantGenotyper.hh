//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Sai Chen <schen6@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <boost/optional.hpp>
#include <vector>

#include "core/Common.hh"
#include "genotyping/SmallVariantGenotype.hh"

namespace ehunter
{

class SmallVariantGenotyper
{
public:
    SmallVariantGenotyper(double haplotypeDepth, AlleleCount expectedAlleleCount)
        : haplotypeDepth_(haplotypeDepth)
        , expectedAlleleCount_((int)expectedAlleleCount)
    {
    }

    /**
     * return the most likely genotype given the read count in reference and alternative allele
     */
    boost::optional<SmallVariantGenotype> genotype(int refCount, int altCount) const;

private:
    /**
     * return a vector of all possible genotypes given the number alleles
     */
    std::vector<SmallVariantGenotype> getPossibleGenotypes(int numAlleles) const;

    /**
     * return genotype likelihood of the given genotype
     * @param currentGenotyper the given genotype
     * @param readCounts Read count vector for each allele (in order)
     */
    double genotypeLikelihood(const SmallVariantGenotype& currentGenotype, int refCount, int altCount) const;

    /**
     * expected depth for one allele
     */
    double haplotypeDepth_;

    /**
     * the expected number of alleles in a genotype (ploidy)
     */
    int expectedAlleleCount_;

    /**
     * hard-coded parameters
     */
    double errorRate_ = 0.05;
};

}
