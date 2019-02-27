//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
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

#include "genotyping/RepeatGenotyper.hh"

#include <algorithm>
#include <cassert>
#include <iostream>

#include <boost/math/distributions/poisson.hpp>

#include "genotyping/RepeatLength.hh"
#include "genotyping/ShortRepeatGenotyper.hh"

namespace ehunter
{

using boost::optional;
using std::map;
using std::string;
using std::vector;

static CountTable combineFlankingAndInrepeatReads(
    int maxNumUnitsInRead, const CountTable& flankingCounts, const CountTable& inrepeatCounts)
{
    int maxNumReadsToTransfer = 5;
    CountTable updatedFlankingCounts = flankingCounts;
    for (int numUnits = maxNumUnitsInRead; numUnits != 0; --numUnits)
    {
        const int count = inrepeatCounts.countOf(numUnits);
        const int countToTransfer = std::min(count, maxNumReadsToTransfer);

        for (int counter = 0; counter != countToTransfer; ++counter)
        {
            updatedFlankingCounts.incrementCountOf(numUnits);
        }

        maxNumReadsToTransfer -= countToTransfer;
        if (maxNumReadsToTransfer == 0)
        {
            break;
        }
    }

    return updatedFlankingCounts;
}

optional<RepeatGenotype> RepeatGenotyper::genotypeRepeat(const vector<int32_t>& alleleSizeCandidates) const
{
    if (alleleSizeCandidates.empty())
    {
        return optional<RepeatGenotype>();
    }

    const CountTable countsOfFlankingReadsForShortRepeatGenotyper
        = combineFlankingAndInrepeatReads(maxNumUnitsInRead_, countsOfFlankingReads_, countsOfInrepeatReads_);

    ShortRepeatGenotyper shortRepeatGenotyper(repeatUnitLen_, maxNumUnitsInRead_, propCorrectMolecules_);

    const int repeatReadCount
        = ::ehunter::countFullLengthRepeatReads(maxNumUnitsInRead_, countsOfFlankingReads_, countsOfInrepeatReads_);

    if (expectedAlleleCount_ == AlleleCount::kOne)
    {
        RepeatGenotype genotype = shortRepeatGenotyper.genotypeRepeatWithOneAllele(
            countsOfFlankingReadsForShortRepeatGenotyper, countsOfSpanningReads_, alleleSizeCandidates);

        const bool isSpanningAllele = countsOfSpanningReads_.countOf(genotype.longAlleleSizeInUnits()) != 0;

        if (!isSpanningAllele && repeatReadCount != 0)
        {
            extendGenotypeWhenOneAlleleIsRepeat(genotype, repeatReadCount);
        }
        else if (!isSpanningAllele)
        {
            extendGenotypeWhenOneAlleleIsFlanking(genotype);
        }
        else
        {
            assert(countsOfSpanningReads_.countOf(genotype.longAlleleSizeInUnits()));
        }

        return genotype;
    }

    assert(expectedAlleleCount_ == AlleleCount::kTwo);

    RepeatGenotype genotype = shortRepeatGenotyper.genotypeRepeatWithTwoAlleles(
        countsOfFlankingReadsForShortRepeatGenotyper, countsOfSpanningReads_, alleleSizeCandidates);

    const bool shortAlleleIsSpanning = countsOfSpanningReads_.countOf(genotype.shortAlleleSizeInUnits()) != 0;
    const bool longAlleleIsSpanning = countsOfSpanningReads_.countOf(genotype.longAlleleSizeInUnits()) != 0;
    // const int repeatReadCount = countFullLengthRepeatReads();

    if (!longAlleleIsSpanning && !shortAlleleIsSpanning && repeatReadCount != 0)
    {
        extendGenotypeWhenBothAllelesAreRepeat(genotype, repeatReadCount);
    }
    else if (!longAlleleIsSpanning && repeatReadCount != 0)
    {
        extendGenotypeWhenOneAlleleIsRepeat(genotype, repeatReadCount);
    }
    else if (shortAlleleIsSpanning && longAlleleIsSpanning)
    {
        // Nothing needs to be done.
    }
    else if (shortAlleleIsSpanning)
    {
        // assert(countsOfFlankingReads_.countOf(genotype.longAlleleSizeInUnits()));
        extendGenotypeWhenOneAlleleIsFlanking(genotype);
    }
    else
    {
        // Both alleles must be flanking.
        // assert(countsOfFlankingReads_.countOf(genotype.shortAlleleSizeInUnits()));
        // assert(countsOfFlankingReads_.countOf(genotype.longAlleleSizeInUnits()));

        extendGenotypeWhenBothAllelesAreFlanking(genotype);
    }

    return genotype;
}

void RepeatGenotyper::extendGenotypeWhenBothAllelesAreFlanking(RepeatGenotype& genotype) const
{
    int32_t flankingAlleleSize, flankingAlleleCiLower, flankingAlleleCiUpper;
    estimateFlankingAlleleSize(flankingAlleleSize, flankingAlleleCiLower, flankingAlleleCiUpper);

    // genotype.setLongAlleleSizeInUnits(flankingAlleleSize);
    flankingAlleleCiLower = std::min(genotype.shortAlleleSizeInUnits(), flankingAlleleCiLower);
    flankingAlleleCiUpper = std::max(genotype.longAlleleSizeInUnits(), flankingAlleleCiUpper);

    genotype.setLongAlleleSizeInUnitsCi(flankingAlleleCiLower, flankingAlleleCiUpper);

    // genotype.setShortAlleleSizeInUnits(flankingAlleleSize);
    genotype.setShortAlleleSizeInUnitsCi(flankingAlleleCiLower, flankingAlleleCiUpper);
}

void RepeatGenotyper::extendGenotypeWhenOneAlleleIsFlanking(RepeatGenotype& genotype) const
{
    int32_t flankingAlleleSize, flankingAlleleCiLower, flankingAlleleCiUpper;
    estimateFlankingAlleleSize(flankingAlleleSize, flankingAlleleCiLower, flankingAlleleCiUpper);

    // genotype.setLongAlleleSizeInUnits(flankingAlleleSize);
    flankingAlleleCiLower = std::min(genotype.shortAlleleSizeInUnits(), flankingAlleleCiLower);
    flankingAlleleCiUpper = std::max(genotype.longAlleleSizeInUnits(), flankingAlleleCiUpper);

    genotype.setLongAlleleSizeInUnitsCi(flankingAlleleCiLower, flankingAlleleCiUpper);
}

void RepeatGenotyper::extendGenotypeWhenOneAlleleIsRepeat(RepeatGenotype& genotype, int numRepeatReads) const
{
    assert(numRepeatReads);

    int32_t longAlleleSize, longAlleleSizeLowerBound, longAlleleSizeUpperBound;
    estimateRepeatAlleleSize(numRepeatReads, longAlleleSize, longAlleleSizeLowerBound, longAlleleSizeUpperBound);
    genotype.setLongAlleleSizeInUnits(longAlleleSize);
    genotype.setLongAlleleSizeInUnitsCi(longAlleleSizeLowerBound, longAlleleSizeUpperBound);
}

void RepeatGenotyper::extendGenotypeWhenBothAllelesAreRepeat(RepeatGenotype& genotype, int numRepeatReads) const
{
    assert(numRepeatReads);

    // Calculate CI for the long allele.
    int32_t longAlleleSize, longAlleleSizeLowerBound, longAlleleSizeUpperBound;
    estimateRepeatAlleleSize(numRepeatReads, longAlleleSize, longAlleleSizeLowerBound, longAlleleSizeUpperBound);
    genotype.setLongAlleleSizeInUnits(longAlleleSize);
    genotype.setLongAlleleSizeInUnitsCi(longAlleleSizeLowerBound, longAlleleSizeUpperBound);

    // Calculate CI for the short allele.
    int32_t shortAlleleSize, shortAlleleSizeLowerBound, shortAlleleSizeUpperBound;
    estimateRepeatAlleleSize(numRepeatReads / 2, shortAlleleSize, shortAlleleSizeLowerBound, shortAlleleSizeUpperBound);
    genotype.setShortAlleleSizeInUnits(shortAlleleSize);
    genotype.setShortAlleleSizeInUnitsCi(shortAlleleSizeLowerBound, shortAlleleSizeUpperBound);
}

void RepeatGenotyper::estimateRepeatAlleleSize(
    int32_t numIrrs, int32_t& size, int32_t& sizeLowerBound, int32_t& sizeUpperBound) const
{
    const int32_t readLength = repeatUnitLen_ * maxNumUnitsInRead_;
    estimateRepeatLen(numIrrs, readLength, haplotypeDepth_, size, sizeLowerBound, sizeUpperBound);

    size /= repeatUnitLen_;
    sizeLowerBound /= repeatUnitLen_;
    sizeUpperBound /= repeatUnitLen_;
}

void RepeatGenotyper::estimateFlankingAlleleSize(
    int32_t& flankingAlleleSize, int32_t& flankingAlleleCiLower, int32_t& flankingAlleleCiUpper) const
{
    const int32_t readLength = repeatUnitLen_ * maxNumUnitsInRead_;

    const vector<int32_t>& spanningSizes = countsOfSpanningReads_.getElementsWithNonzeroCounts();
    const int32_t longestSpanning
        = spanningSizes.empty() ? 0 : *std::max_element(spanningSizes.begin(), spanningSizes.end());

    int32_t numFlankingReadsLongerThanSpanning = 0;
    for (int32_t repeatSize : countsOfFlankingReads_.getElementsWithNonzeroCounts())
    {
        if (repeatSize > longestSpanning)
        {
            numFlankingReadsLongerThanSpanning += countsOfFlankingReads_.countOf(repeatSize);
        }
    }

    // Haplotype depth should be twice as high because flanking reads come from both flanks of the repeat.
    estimateRepeatLen(
        numFlankingReadsLongerThanSpanning, readLength, 2 * haplotypeDepth_, flankingAlleleSize, flankingAlleleCiLower,
        flankingAlleleCiUpper);

    // estimateRepeatLen adds read length to size estimates so we need to subtract it out.
    flankingAlleleSize -= readLength;
    flankingAlleleCiLower -= readLength;
    flankingAlleleCiUpper -= readLength;

    flankingAlleleSize = flankingAlleleSize / repeatUnitLen_ + longestSpanning + 1;
    flankingAlleleCiLower = flankingAlleleCiLower / repeatUnitLen_ + longestSpanning + 1;
    flankingAlleleCiUpper = flankingAlleleCiUpper / repeatUnitLen_ + longestSpanning + 1;

    // Repeat must be at least at long as the longest flanking read.
    const vector<int32_t>& flankingSizes = countsOfFlankingReads_.getElementsWithNonzeroCounts();
    const int32_t longestFlanking = *std::max_element(flankingSizes.begin(), flankingSizes.end());

    flankingAlleleCiLower = std::max(flankingAlleleCiLower, longestFlanking);
    flankingAlleleSize = std::max(flankingAlleleSize, longestFlanking);
    flankingAlleleCiUpper = std::max(flankingAlleleCiUpper, longestFlanking);

    // Repeat estimated from flanking reads cannot be longer than the read length.
    flankingAlleleCiLower = std::min(flankingAlleleCiLower, maxNumUnitsInRead_);
    flankingAlleleSize = std::min(flankingAlleleSize, maxNumUnitsInRead_);
    flankingAlleleCiUpper = std::min(flankingAlleleCiUpper, maxNumUnitsInRead_);
}

int RepeatGenotyper::countFullLengthRepeatReads() const
{
    const double kMinProportionAlignedBases = 0.93;
    const int fullLengthSizeCutoff = static_cast<int>(std::round(maxNumUnitsInRead_ * kMinProportionAlignedBases));

    int repeatReadCount = 0;
    for (int numUnits : countsOfInrepeatReads_.getElementsWithNonzeroCounts())
    {
        if (numUnits >= fullLengthSizeCutoff)
        {
            repeatReadCount += countsOfInrepeatReads_.countOf(numUnits);
        }
    }

    return repeatReadCount;
}

static int depthBasedCountOfInrepeatReads(
    int maxNumUnitsInRead, const CountTable& countsOfFlankingReads, const CountTable& countsOfInrepeatReads)
{
    const int kNumFlanks = 2;
    const double kPropLowConfidenceFlank = 0.1;
    const double kPropHighConfidenceFlank = 1.0 - kPropLowConfidenceFlank;

    const int maxUnitsFromLowConfidenceFlank
        = static_cast<int>(std::round(maxNumUnitsInRead * kPropHighConfidenceFlank));

    int numPutativeIrrs = 0;
    for (const auto& numUnitsSpannedAndCount : countsOfInrepeatReads)
    {
        int numUnitsSpanned = numUnitsSpannedAndCount.first;
        int readCount = numUnitsSpannedAndCount.second;

        if (numUnitsSpanned >= maxUnitsFromLowConfidenceFlank)
        {
            numPutativeIrrs += readCount;
        }
    }

    int numFlankingReads = 0;
    for (const auto& numUnitsSpannedAndCount : countsOfFlankingReads)
    {
        int numUnitsSpanned = numUnitsSpannedAndCount.first;
        int readCount = numUnitsSpannedAndCount.second;

        if (numUnitsSpanned < maxUnitsFromLowConfidenceFlank)
        {
            numFlankingReads += readCount;
        }
    }
    if (numFlankingReads == 0) return 0;

    const double estimatedDepth = numFlankingReads / (kNumFlanks * kPropHighConfidenceFlank);
    const double expectedNumLowconfidenceFlankingReads = kNumFlanks * kPropLowConfidenceFlank * estimatedDepth;

    boost::math::poisson_distribution<> lowconfidenceFlankingDistro(expectedNumLowconfidenceFlankingReads);
    const double probability = boost::math::cdf(lowconfidenceFlankingDistro, numPutativeIrrs);

    const double kProbabilityCutoff = 0.95;
    return (probability >= kProbabilityCutoff ? numPutativeIrrs : 0);
}

static int lengthBasedCountOfInrepeatReads(int maxNumUnitsInRead, const CountTable& countsOfInrepeatReads)
{
    const double kPropForFullLength = 0.96;
    const int minNumUnitsForFullLength = static_cast<int>(std::round(maxNumUnitsInRead * kPropForFullLength));

    int numIrrs = 0;
    for (const auto& numUnitsSpannedAndCount : countsOfInrepeatReads)
    {
        int numUnitsSpanned = numUnitsSpannedAndCount.first;
        int readCount = numUnitsSpannedAndCount.second;

        if (numUnitsSpanned >= minNumUnitsForFullLength)
        {
            numIrrs += readCount;
        }
    }

    return numIrrs;
}

int countFullLengthRepeatReads(
    int maxNumUnitsInRead, const CountTable& countsOfFlankingReads, const CountTable& countsOfInrepeatReads)
{
    const int lengthBasedCount = lengthBasedCountOfInrepeatReads(maxNumUnitsInRead, countsOfInrepeatReads);
    const int depthBasedCount = depthBasedCountOfInrepeatReads(maxNumUnitsInRead, countsOfFlankingReads, countsOfInrepeatReads);
    return std::max(lengthBasedCount, depthBasedCount);
}

}
