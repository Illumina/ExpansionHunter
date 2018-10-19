//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "genotyping/RepeatGenotyper.hh"

#include <algorithm>
#include <cassert>
#include <iostream>

#include "genotyping/RepeatLength.hh"
#include "genotyping/ShortRepeatGenotyper.hh"

using boost::optional;
using std::cerr;
using std::endl;
using std::map;
using std::string;
using std::vector;

optional<RepeatGenotype> RepeatGenotyper::genotypeRepeat(const vector<int32_t>& alleleSizeCandidates) const
{
    if (alleleSizeCandidates.empty())
    {
        return optional<RepeatGenotype>();
    }

    ShortRepeatGenotyper shortRepeatGenotyper(repeatUnitLen_, maxNumUnitsInRead_, propCorrectMolecules_);

    if (expectedAlleleCount_ == AlleleCount::kOne)
    {
        RepeatGenotype genotype = shortRepeatGenotyper.genotypeRepeatWithOneAllele(
            countsOfFlankingReads_, countsOfSpanningReads_, alleleSizeCandidates);

        const bool isSpanningAllele = countsOfSpanningReads_.countOf(genotype.longAlleleSizeInUnits()) != 0;
        const int repeatReadCount = countFullLengthRepeatReads();

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
        countsOfFlankingReads_, countsOfSpanningReads_, alleleSizeCandidates);

    const bool shortAlleleIsSpanning = countsOfSpanningReads_.countOf(genotype.shortAlleleSizeInUnits()) != 0;
    const bool longAlleleIsSpanning = countsOfSpanningReads_.countOf(genotype.longAlleleSizeInUnits()) != 0;
    const int repeatReadCount = countFullLengthRepeatReads();

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
        assert(countsOfFlankingReads_.countOf(genotype.longAlleleSizeInUnits()));
        extendGenotypeWhenOneAlleleIsFlanking(genotype);
    }
    else
    {
        // Both alleles must be flanking.
        assert(countsOfFlankingReads_.countOf(genotype.shortAlleleSizeInUnits()));
        assert(countsOfFlankingReads_.countOf(genotype.longAlleleSizeInUnits()));

        extendGenotypeWhenBothAllelesAreFlanking(genotype);
    }

    return genotype;
}

void RepeatGenotyper::extendGenotypeWhenBothAllelesAreFlanking(RepeatGenotype& genotype) const
{
    int32_t flankingAlleleSize, flankingAlleleCiLower, flankingAlleleCiUpper;
    estimateFlankingAlleleSize(flankingAlleleSize, flankingAlleleCiLower, flankingAlleleCiUpper);

    genotype.setLongAlleleSizeInUnits(flankingAlleleSize);
    genotype.setLongAlleleSizeInUnitsCi(flankingAlleleCiLower, flankingAlleleCiUpper);

    genotype.setShortAlleleSizeInUnits(flankingAlleleSize);
    genotype.setShortAlleleSizeInUnitsCi(flankingAlleleCiLower, flankingAlleleCiUpper);
}

void RepeatGenotyper::extendGenotypeWhenOneAlleleIsFlanking(RepeatGenotype& genotype) const
{
    int32_t flankingAlleleSize, flankingAlleleCiLower, flankingAlleleCiUpper;
    estimateFlankingAlleleSize(flankingAlleleSize, flankingAlleleCiLower, flankingAlleleCiUpper);

    genotype.setLongAlleleSizeInUnits(flankingAlleleSize);
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
    const int fullLengthSizeCutoff = maxNumUnitsInRead_ - 2;

    int repeatReadCount = 0;
    for (int numUnits : countsOfRepeatReads_.getElementsWithNonzeroCounts())
    {
        if (numUnits >= fullLengthSizeCutoff)
        {
            repeatReadCount += countsOfRepeatReads_.countOf(numUnits);
        }
    }

    return repeatReadCount;
}
