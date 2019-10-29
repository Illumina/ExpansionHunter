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

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "common/Common.hh"
#include "common/CountTable.hh"
#include "common/Parameters.hh"
#include "genotyping/RepeatGenotype.hh"
#include "locus_spec/LocusSpec.hh"

namespace ehunter
{

class RepeatGenotyper
{
public:
    RepeatGenotyper(
        double haplotypeDepth, AlleleCount expectedAlleleCount, int32_t repeatUnitLen, int32_t maxNumUnitsInRead,
        double propCorrectMolecules, const CountTable& countsOfSpanningReads, const CountTable& countsOfFlankingReads,
        const CountTable& countsOfRepeatReads, int countOfInrepeatReadPairs)
        : expectedAlleleCount_(expectedAlleleCount)
        , repeatUnitLen_(repeatUnitLen)
        , maxNumUnitsInRead_(maxNumUnitsInRead)
        , haplotypeDepth_(haplotypeDepth)
        , propCorrectMolecules_(propCorrectMolecules)
        , countsOfSpanningReads_(countsOfSpanningReads)
        , countsOfFlankingReads_(countsOfFlankingReads)
        , countsOfInrepeatReads_(countsOfRepeatReads)
        , countOfInrepeatReadPairs_(countOfInrepeatReadPairs)
    {
    }

    boost::optional<RepeatGenotype> genotypeRepeat(const std::vector<int32_t>& alleleSizeCandidates) const;

    // The methods below are exposed for unit-testing purposes.
    // When both alleles are longer than the read length we cannot, in general, know which allele a given in-repeat read
    // originated from. To account for this we compute confidence intervals corresponding to two extreme cases of
    // partitioning IRR proportions between the two alleles (0.5/0.5 and 0/1.0) and compute widest possible confidence
    // intervals based on these estimates.
    void extendGenotypeWhenBothAllelesAreRepeat(RepeatGenotype& genotype, int numRepeatReads) const;
    void extendGenotypeWhenOneAlleleIsRepeat(RepeatGenotype& genotype, int numRepeatReads) const;
    void extendGenotypeWhenBothAllelesAreFlanking(RepeatGenotype& genotype) const;
    void extendGenotypeWhenOneAlleleIsFlanking(RepeatGenotype& genotype) const;

private:
    int calculateLongestSpanning() const;
    int countFlankingReadsLongerThanSpanning() const;

    void estimateFlankingAlleleSize(
        int32_t& flankingAlleleSize, int32_t& flankingAlleleCiLower, int32_t& flankingAlleleCiUpper) const;
    void estimateRepeatAlleleSize(
        int32_t numIrrs, int32_t& flankingAlleleSize, int32_t& flankingAlleleCiLower,
        int32_t& flankingAlleleCiUpper) const;

    int countFullLengthRepeatReads() const;

    const AlleleCount expectedAlleleCount_;
    const int32_t repeatUnitLen_;
    const int32_t maxNumUnitsInRead_;
    const double haplotypeDepth_;
    const double propCorrectMolecules_;
    const CountTable countsOfSpanningReads_;
    const CountTable countsOfFlankingReads_;
    const CountTable countsOfInrepeatReads_;
    int countOfInrepeatReadPairs_;
};

int countFullLengthRepeatReads(
    int maxNumUnitsInRead, const CountTable& countsOfFlankingReads, const CountTable& countsOfInrepeatReads);

}
