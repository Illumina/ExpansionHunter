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
#include "region_spec/LocusSpecification.hh"

namespace ehunter
{

class RepeatGenotyper
{
public:
    RepeatGenotyper(
        double haplotypeDepth, AlleleCount expectedAlleleCount, int32_t repeatUnitLen, int32_t maxNumUnitsInRead,
        double propCorrectMolecules, const CountTable& countsOfSpanningReads, const CountTable& countsOfFlankingReads,
        const CountTable& countsOfRepeatReads)
        : expectedAlleleCount_(expectedAlleleCount)
        , repeatUnitLen_(repeatUnitLen)
        , maxNumUnitsInRead_(maxNumUnitsInRead)
        , haplotypeDepth_(haplotypeDepth)
        , propCorrectMolecules_(propCorrectMolecules)
        , countsOfSpanningReads_(countsOfSpanningReads)
        , countsOfFlankingReads_(countsOfFlankingReads)
        , countsOfInrepeatReads_(countsOfRepeatReads)
    {
    }

    boost::optional<RepeatGenotype> genotypeRepeat(const std::vector<int32_t>& alleleSizeCandidates) const;

private:
    void extendGenotypeWhenBothAllelesAreRepeat(RepeatGenotype& genotype, int numRepeatReads) const;
    void extendGenotypeWhenOneAlleleIsRepeat(RepeatGenotype& genotype, int numRepeatReads) const;
    void extendGenotypeWhenBothAllelesAreFlanking(RepeatGenotype& genotype) const;
    void extendGenotypeWhenOneAlleleIsFlanking(RepeatGenotype& genotype) const;

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
};

int countFullLengthRepeatReads(
    int maxNumUnitsInRead, const CountTable& countsOfFlankingReads, const CountTable& countsOfInrepeatReads);

}
