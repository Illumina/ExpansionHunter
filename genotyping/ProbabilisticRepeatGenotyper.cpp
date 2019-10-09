//
// Expansion Hunter
// Copyright (c) 2019 Illumina, Inc.
//
// Author: Konrad Scheffler <kscheffler@illumina.com>
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

#include "genotyping/ProbabilisticRepeatGenotyper.hh"
#include "stats/LogSumUtil.hh"

#include <algorithm>
#include <math.h>

namespace ehunter
{

using std::map;
using std::string;
using std::vector;

const static double logOneHalf(std::log(0.5));

void ProbabilisticGenotypeScoreSet::addGenotypeScore(const RepeatGenotype genotype, double score)
{
    // The design is to insert everything in arbitrary order, sort once after insertion stage,
    // then keep the order unchanged. Keeping track of whether the vector is in sorted state is just a precaution in
    // case a client deviates from this design, so we don't bother implementing this in the most efficient way.
    isSorted_ = false;
    isNormalized_ = false;
    scoredGenotypes_.emplace_back(ScoredGenotype(genotype, score));
}

optional<RepeatGenotype> ProbabilisticGenotypeScoreSet::bestGenotype()
{
    optional<RepeatGenotype> bestGenotype;
    if (!scoredGenotypes_.empty())
    {
        bestGenotype = scoredGenotypes_[0].genotype;
    }
    return bestGenotype;
}

void ProbabilisticGenotypeScoreSet::normalize()
{
    if (scoredGenotypes_.empty() || isNormalized_)
        return;
    if (!isSorted_)
    {
        std::sort(scoredGenotypes_.begin(), scoredGenotypes_.end());
        isSorted_ = true;
    }

    // Calculate normalizing constant: this is the sum over all linear-domain probabilities,
    // calculated by first normalizing by max to alleviate underflow (see logSumUtil):
    auto it = std::begin(scoredGenotypes_);
    double bestScore(it->score);
    ++it;
    double totalProb(0.0);
    for (auto end = std::end(scoredGenotypes_); it != end; ++it)
    {
        totalProb += std::exp(it->score - bestScore);
    }
    double logNormalizingConstant(bestScore + log1p_switch(totalProb));

    // Finally, do the actual normalization:
    for (auto& elem : scoredGenotypes_)
    {
        elem.score -= logNormalizingConstant;
    }
    isNormalized_ = true;
}

void ProbabilisticGenotypeScoreSet::constructCredibleInterval(double credibleIntervalSize)
{
    if (scoredGenotypes_.empty())
        return;
    if (!isNormalized_)
        normalize();

    double cumulativeProb(0.0);
    int32_t shortMin(scoredGenotypes_[0].genotype.shortAlleleSizeInUnits());
    int32_t shortMax(scoredGenotypes_[0].genotype.shortAlleleSizeInUnits());
    int32_t longMin(scoredGenotypes_[0].genotype.longAlleleSizeInUnits());
    int32_t longMax(scoredGenotypes_[0].genotype.longAlleleSizeInUnits());

    // iterate through set in order of descending score, until cumulative score exceeds interval size:
    for (const auto& elem : scoredGenotypes_)
    {
        int32_t shortAllele(elem.genotype.shortAlleleSizeInUnits());
        int32_t longAllele(elem.genotype.shortAlleleSizeInUnits());

        if (shortAllele < shortMin)
            shortMin = shortAllele;
        if (shortAllele > shortMax)
            shortMax = shortAllele;
        if (longAllele < longMin)
            longMin = longAllele;
        if (longAllele > longMax)
            longMax = longAllele;

        cumulativeProb += std::exp(elem.score);
        if (cumulativeProb >= credibleIntervalSize)
            break;
    }

    scoredGenotypes_[0].genotype.setShortAlleleSizeInUnitsCi(shortMin, shortMax);
    scoredGenotypes_[0].genotype.setLongAlleleSizeInUnitsCi(longMin, longMax);
}

vector<vector<double>> ProbabilisticRepeatGenotyper::scoreReadsAgainstAlleles() const
{
    vector<vector<double>> readAlleleScores; // allele score vectors for all reads
    /*
    for each read:
        for each alignment:
            for each allele:
                update readAllele score using alignment score and gap penalty
    */
    for (auto const& read : readSummaries_)
    {
        vector<double> alleleScores(maxAlleleSize_ + 1, -INFINITY); // allele scores for this read
        for (auto alignment : read.alignments())
        {
            // primary allele associated with this alignment:
            int32_t primaryAllele(alignment.numUnits());

            // 1) Get stutter-free likelihood for read given this alignment. To be adjusted for stutter alleles.
            // lnP(read | graphAlignment, allele) :
            double readLnLGivenPrimaryAllele(alignment.score());

            // 2) Get alignment prior given stutter-free allele. To be adjusted for stutter/alternate alleles.
            int32_t regionLength(adjustedRegionSize_ + primaryAllele * repeatUnitLen_);
            int32_t numPossibleAlignmentPositions(alignment.clippedReadLength() + regionLength - 1);
            int32_t numActualAlignmentPositions(1);
            // lnP(graphAlignment | allele) :
            double graphAlignmentLnProbGivenAllele(
                -std::log(numPossibleAlignmentPositions)); // 1/len for non-repeat reads

            // Combine (1) and (2) into primary allele score:
            // lnP(read, graphAlignment | allele) :
            double readGraphAlignmentJointLnProb(readLnLGivenPrimaryAllele + graphAlignmentLnProbGivenAllele);
            alleleScores[primaryAllele] = getLogSum(alleleScores[primaryAllele], readGraphAlignmentJointLnProb);

            // Loop through non-primary alleles:
            double readLnLGivenStutterAllele(readLnLGivenPrimaryAllele);
            int32_t numPossibleShortAlignmentPositions(numPossibleAlignmentPositions);
            int32_t numPossibleLongAlignmentPositions(numPossibleAlignmentPositions);
            // numActualShortAlignmentPositions is always 1 even for in-repeat reads,
            // because the primary alignment is the shortest one consistent with the read
            int32_t numActualLongAlignmentPositions(numActualAlignmentPositions);
            for (int32_t shortAllele = primaryAllele - 1, longAllele = primaryAllele + 1;
                 (shortAllele >= 0) || (longAllele <= maxAlleleSize_); --shortAllele, ++longAllele)
            {
                readLnLGivenStutterAllele += stutterPenalty_;

                if (shortAllele >= 0)
                {
                    numPossibleShortAlignmentPositions -= repeatUnitLen_;
                    graphAlignmentLnProbGivenAllele = -std::log(numPossibleShortAlignmentPositions);
                    alleleScores[shortAllele] = getLogSum(
                        alleleScores[shortAllele], readLnLGivenStutterAllele + graphAlignmentLnProbGivenAllele);
                }

                if (longAllele <= maxAlleleSize_)
                {
                    double readLnLGivenAllele(
                        alignment.isSpanning() ? readLnLGivenStutterAllele : readLnLGivenPrimaryAllele);
                    numPossibleLongAlignmentPositions += repeatUnitLen_;
                    graphAlignmentLnProbGivenAllele = -std::log(numPossibleLongAlignmentPositions);
                    if (alignment.isRepeat())
                    {
                        ++numActualLongAlignmentPositions;
                        graphAlignmentLnProbGivenAllele += std::log(numActualLongAlignmentPositions);
                    }
                    alleleScores[longAllele]
                        = getLogSum(alleleScores[longAllele], readLnLGivenAllele + graphAlignmentLnProbGivenAllele);
                }
            }
        } // for alignment
        readAlleleScores.push_back(alleleScores);
    } // for read
    return readAlleleScores;
}

double ProbabilisticRepeatGenotyper::alleleBias(double alleleOne, double alleleTwo) const
{
    double alleleOneOverlapOpportunities(alleleOne * repeatUnitLen_ + expectedReadLength_ - 1);
    double alleleTwoOverlapOpportunities(alleleTwo * repeatUnitLen_ + expectedReadLength_ - 1);
    return alleleOneOverlapOpportunities / (alleleOneOverlapOpportunities + alleleTwoOverlapOpportunities);
}

ProbabilisticGenotypeScoreSet
ProbabilisticRepeatGenotyper::scoreGenotypes(const vector<vector<double>>& readAlleleScores) const
{
    ProbabilisticGenotypeScoreSet lnPosteriors;
    if (ploidy_ == AlleleCount::kOne) // haploid case
    {
        /*
         *  for each haploid genotype:
         *      calculate unnormalized genotype posterior, keeping track of best score
         */
        for (int32_t allele = 0; allele <= maxAlleleSize_; ++allele)
        {
            RepeatGenotype genotype(repeatUnitLen_, { allele });
            double lnPrior(
                0.0); // A constant value here gives a uniform prior TODO: implement a more sophisticated prior
            double lnL(0.0);
            for (int32_t readIdx = 0; readIdx < (int32_t)readAlleleScores.size(); ++readIdx)
            {
                double lnLGivenMismap(readSummaries_[readIdx].readLength() * randomBasePenalty_);
                lnL += getLogSum(
                    lnLGivenMismap + mismapLnPrior_, readAlleleScores[readIdx][allele] + correctmapLnPrior_);
            }
            lnPosteriors.addGenotypeScore(genotype, lnPrior + lnL);
        }
    } // end haploid case
    else if (ploidy_ == AlleleCount::kTwo) // diploid case
    {
        /*
         *  for each diploid genotype:
         *      calculate unnormalized genotype posterior, keeping track of best score
         */
        for (int32_t alleleOne = 0; alleleOne <= maxAlleleSize_; ++alleleOne)
        {
            for (int32_t alleleTwo = alleleOne; alleleTwo <= maxAlleleSize_; ++alleleTwo)
            {
                RepeatGenotype genotype(repeatUnitLen_, { alleleOne, alleleTwo });
                double lnPrior(
                    alleleOne == alleleTwo
                        ? logOneHalf
                        : 0.0); // Make homozygotes half as likely as heterozygotes; prior is uniform in other respects
                // TODO: implement a more sophisticated prior
                double lnL(0.0);
                for (int32_t readIdx = 0; readIdx < (int32_t)readAlleleScores.size(); ++readIdx)
                {
                    double lnLGivenMismap(readSummaries_[readIdx].readLength() * randomBasePenalty_);
                    double alleleOneSampleProb(alleleBias(alleleOne, alleleTwo));
                    double lnLGivenCorrectlyMapped = getLogSum(
                        readAlleleScores[readIdx][alleleOne] + std::log(alleleOneSampleProb),
                        readAlleleScores[readIdx][alleleTwo] + std::log(1.0 - alleleOneSampleProb));
                    lnL += getLogSum(lnLGivenMismap + mismapLnPrior_, lnLGivenCorrectlyMapped + correctmapLnPrior_);
                }
                lnPosteriors.addGenotypeScore(genotype, lnPrior + lnL);
            }
        }
    } // end diploid case

    return lnPosteriors;
}

optional<RepeatGenotype> ProbabilisticRepeatGenotyper::genotypeRepeat(double credibleIntervalSize) const
{
    assert(credibleIntervalSize > 0.0 && credibleIntervalSize < 1.0);

    auto readScores = scoreReadsAgainstAlleles();
    auto genotypeScores = scoreGenotypes(readScores);
    genotypeScores.constructCredibleInterval(credibleIntervalSize);
    return genotypeScores.bestGenotype();
}
}
