//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Sai Chen <schen6@illumina.com>
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

#include "genotyping/SmallVariantGenotyper.hh"
#include <boost/math/distributions/poisson.hpp>
#include <numeric>

using boost::optional;
using boost::math::poisson_distribution;
using std::vector;

namespace ehunter
{

boost::optional<SmallVariantGenotype> SmallVariantGenotyper::genotype(int refCount, int altCount) const
{
    if (expectedAlleleCount_ > 2)
    {
        throw std::runtime_error(
            "Ploidy = " + std::to_string(expectedAlleleCount_) + ", but only ploidy <=2 is supported at this time.");
    }

    if (refCount < 0 || altCount < 0)
    {
        throw std::runtime_error("Invalid read counts: " + std::to_string(refCount) + " " + std::to_string(altCount));
    }

    auto possibleGenotypes = getPossibleGenotypes(expectedAlleleCount_);

    const int totalReadCount = refCount + altCount;
    if (totalReadCount == 0) // missing genotype
    {
        return {};
    }

    const unsigned genotypeCount(possibleGenotypes.size());
    if (genotypeCount == 0)
    {
        return {};
    }

    unsigned mostLikelyGenotypeIndex(0);
    double bestLikelihood(0);
    for (unsigned genotypeIndex(0); genotypeIndex < genotypeCount; ++genotypeIndex)
    {
        const double currentLikelihood = genotypeLikelihood(possibleGenotypes[genotypeIndex], refCount, altCount);
        if ((genotypeIndex == 0) or (currentLikelihood > bestLikelihood))
        {
            bestLikelihood = currentLikelihood;
            mostLikelyGenotypeIndex = genotypeIndex;
        }
    }

    return possibleGenotypes[mostLikelyGenotypeIndex];
}

vector<SmallVariantGenotype> SmallVariantGenotyper::getPossibleGenotypes(int numAlleles) const
{
    if (numAlleles == 1)
    {
        return { AlleleType::kRef, AlleleType::kAlt };
    }
    else if (numAlleles == 2)
    {
        return { { AlleleType::kRef, AlleleType::kRef },
                 { AlleleType::kRef, AlleleType::kAlt },
                 { AlleleType::kAlt, AlleleType::kAlt } };
    }
    else
    {
        throw std::runtime_error("The given number of alleles is not supported: " + std::to_string(numAlleles));
    }
}

double
SmallVariantGenotyper::genotypeLikelihood(const SmallVariantGenotype& currentGenotype, int refCount, int altCount) const
{
    const poisson_distribution<> errorDistribution(errorRate_);

    const bool isHomozygous = currentGenotype.isHomRef() || currentGenotype.isHomAlt();
    const int copyNumberOfExistingAllele = isHomozygous ? 2 : 1;
    const poisson_distribution<> countDistribution(copyNumberOfExistingAllele * haplotypeDepth_);

    double genotypeLikelihood
        = currentGenotype.isHomRef() ? log(pdf(errorDistribution, altCount)) : log(pdf(countDistribution, altCount));

    genotypeLikelihood
        += currentGenotype.isHomAlt() ? log(pdf(errorDistribution, refCount)) : log(pdf(countDistribution, refCount));

    if (std::isinf(genotypeLikelihood))
    {
        genotypeLikelihood = -std::numeric_limits<double>::max();
    }

    return genotypeLikelihood;
}

}
