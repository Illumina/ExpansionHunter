//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Sai Chen <schen6@illumina.com>,
//         Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#include <boost/optional.hpp>
#include <vector>

#include "common/Common.hh"
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
