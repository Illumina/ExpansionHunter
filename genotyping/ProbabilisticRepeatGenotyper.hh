//
// Expansion Hunter
// Copyright (c) 2019 Illumina, Inc.
//
// Author: Konrad Scheffler <kscheffler@illumina.com>,
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

// Probability-related conventions:
//
// We work in the log-domain where possible, using natural logs (ln). Log-domain quantities are indicated by
// including "ln" in the variable name. When this is absent it means we are working with linear-domain probabilities.
// If (shudder) we ever represent something in base-10 logs, that is indicated by including "log10" in the variable
// name. To avoid potential ambiguity, we never use "log" in variable names. Log-likelihoods are named lnL. These always
// denote the probability of (observed) data given something else. Priors (or log-priors) are named prior (or lnPrior).
// These are always independent of data. Posteriors (or log-posteriors) are named posterior (or lnPosterior). These
// always denote the probability of something given data. We follow the same naming scheme and refer to quantities as
// probabilities even if they are not normalized.

#pragma once

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "classification/AlignmentSummary.hh"
#include "genotyping/RepeatGenotype.hh"

namespace ehunter
{
using boost::optional;
using std::vector;

// Collection of genotypes and associated scores
// TODO: Refactor RepeatGenotype so that scores can be added there directly. Then we can remove this class
//  (member functions can become private members of ProbabilisticRepeatGenotyper and scoredGenotypes_ can become
//  a local variable in ProbabilisticRepeatGenotyper::genotypeRepeat())
class ProbabilisticGenotypeScoreSet
{
public:
    ProbabilisticGenotypeScoreSet()
        : isSorted_(false)
        , isNormalized_(false)
    {
    }

    // Adds a scored genotype to the set.
    void addGenotypeScore(RepeatGenotype genotype, double score);

    // Returns the best-scoring genotype (if the set is non-empty).
    optional<RepeatGenotype> bestGenotype();

    // Normalizes the list of genotype scores.
    void normalize();

    // Constructs credible intervals around every allele in the best genotype.
    void constructCredibleInterval(double credibleIntervalSize);

private:
    struct ScoredGenotype
    {
        ScoredGenotype(RepeatGenotype gt, double s)
            : genotype(std::move(gt))
            , score(s)
        {
        }
        // Overload less-than for sorting in descending order of score:
        bool operator<(ScoredGenotype const& other) const { return score > other.score; }
        RepeatGenotype genotype;
        double score;
    };

    std::vector<ScoredGenotype> scoredGenotypes_;
    bool isSorted_;
    bool isNormalized_;
};

class ProbabilisticRepeatGenotyper
{
public:
    /**
     * @param adjustedRegionSize The expected size (length) of the region covered by the graph,
     * with the STR of interest excised (average over any nuisance STRs).
     * @param expectedReadLength The average read length over the run.
     * @param maxAlleleSize The maximum allele size (number of repeat units) to consider during genotyping.
     * @param stutterPenalty We use a simplistic model in which stutter errors are penalized per repeat unit.
     * The penalty must be a properly normalized log-probability, i.e. a negative number. It will be interpreted as
     * the log-probability of a single repeat unit stutter error and added to the (unnormalized)
     * stutter-free log-probabilities.
     * @param randomBasePenalty The unnormalized log-probability per base in a mismapped read. This should be the same
     * as the score for a clipped base in the alignment scoring scheme. Like stutterPenalty, more negative values denote
     * larger penalties. Unlike stutterPenalty the value is unnormalized and does not need to be negative (so we cannot
     * enforce an assertion - be careful!) For now, use a value of 0. TODO: should we use (5+log(1/4)) instead?
     * @param mismapProb The probability with which reads that have been mapped to this region do not belong here.
     * @param readSummaries A vector of ReadSummary objects (one per read), each containing a vector of alignments.
     * This is the main data structure to be processed.
     */
    ProbabilisticRepeatGenotyper(
        AlleleCount ploidy, int32_t repeatUnitLen, int32_t adjustedRegionSize, int32_t expectedReadLength,
        int32_t maxAlleleSize, double stutterPenalty, double randomBasePenalty, double mismapProb,
        std::vector<ReadSummaryForStr> readSummaries)
        : ploidy_(ploidy)
        , repeatUnitLen_(repeatUnitLen)
        , adjustedRegionSize_(adjustedRegionSize)
        , expectedReadLength_(expectedReadLength)
        , maxAlleleSize_(maxAlleleSize)
        , stutterPenalty_(stutterPenalty)
        , randomBasePenalty_(randomBasePenalty)
        , mismapLnPrior_(std::log(mismapProb))
        , correctmapLnPrior_(std::log(1.0 - mismapProb))
        , readSummaries_(std::move(readSummaries))
    {
        assert(stutterPenalty < 0);
    }

    /**
     * @param credibleIntervalSize e.g. 0.95 for a 95% credible interval
     */
    optional<RepeatGenotype> genotypeRepeat(double credibleIntervalSize = 0.95) const;

private:
    // Scores every read against every allele (and against mismap).
    vector<vector<double>> scoreReadsAgainstAlleles() const;

    // Allele bias calculation: the probability of getting an alleleOne read
    // when sampling from a diploid (alleleOne, alleleTwo) mixture
    double alleleBias(double alleleOne, double alleleTwo) const;

    // Generates and scores a complete list of genotypes. Returns the best-scoring genotype.
    ProbabilisticGenotypeScoreSet scoreGenotypes(const vector<vector<double>>& readAlleleScores) const;

    const AlleleCount ploidy_;
    const int32_t repeatUnitLen_;
    const int32_t adjustedRegionSize_;
    const int32_t expectedReadLength_;
    const int32_t maxAlleleSize_;
    const double stutterPenalty_;
    const double randomBasePenalty_;
    const double mismapLnPrior_;
    const double correctmapLnPrior_;
    const std::vector<ReadSummaryForStr> readSummaries_;
};
}
