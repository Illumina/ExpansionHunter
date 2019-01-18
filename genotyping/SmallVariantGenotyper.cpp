//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Sai Chen <schen6@illumina.com>
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
        return optional<SmallVariantGenotype>();
    }

    optional<SmallVariantGenotype> mostLikelyGenotype;
    double bestLikelihood = -std::numeric_limits<double>::max();
    for (const auto& currentGenotype : possibleGenotypes)
    {
        const double currentLikelihood = genotypeLikelihood(currentGenotype, refCount, altCount);
        if (currentLikelihood > bestLikelihood)
        {
            bestLikelihood = currentLikelihood;
            mostLikelyGenotype = currentGenotype;
        }
    }

    return mostLikelyGenotype;
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
