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

// Defines classes and methods for GenotypeRepeat and haplotype likelihood
// calculations.

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "common/Common.hh"
#include "common/CountTable.hh"
#include "genotyping/RepeatGenotype.hh"

namespace ehunter
{

class QuantifierOfMoleculesGeneratedByAllele
{
public:
    QuantifierOfMoleculesGeneratedByAllele(
        int32_t allele_size_in_units, int32_t max_repeat_size_in_units, double prop_correct_molecules);
    double propMoleculesOfGivenSize(int32_t size_in_units) const;
    double propMoleculesShorterThan(int32_t size_upper_bound_in_units) const;
    double propMoleculesAtLeast(int32_t size_lower_bound_in_units) const;
    int32_t alleleSizeInUnits() const { return allele_size_in_units_; }

private:
    int32_t allele_size_in_units_;
    int32_t max_repeat_size_in_units_;
    double prop_correct_molecules_;
    double norm_factor_;
    int max_deviation_;
};

class ShortRepeatGenotypeLikelihoodEstimator
{
public:
    ShortRepeatGenotypeLikelihoodEstimator(
        int32_t max_repeat_size_in_units, double prop_correct_molecules, std::vector<int32_t> allele_sizes_in_units)
        : max_repeat_size_in_units_(max_repeat_size_in_units)
    {
        for (int32_t allele_size_in_units : allele_sizes_in_units)
        {
            allele_quantifiers_.emplace_back(allele_size_in_units, max_repeat_size_in_units, prop_correct_molecules);
        }
    }

    double CalcFlankingLoglik(int32_t repeat_size_in_units) const;
    double CalcSpanningLoglik(int32_t repeat_size_in_units) const;
    double CalcLogLik(const CountTable& counts_of_flanking_reads, const CountTable& counts_of_spanning_reads) const;

private:
    int max_repeat_size_in_units_;
    std::vector<QuantifierOfMoleculesGeneratedByAllele> allele_quantifiers_;
};

class ShortRepeatGenotyper
{
public:
    ShortRepeatGenotyper(int32_t repeatUnitLen, int32_t max_repeat_size_in_units, double prop_correct_molecules)
        : repeatUnitLen_(repeatUnitLen)
        , max_repeat_size_in_units_(max_repeat_size_in_units)
        , prop_correct_molecules_(prop_correct_molecules)
    {
    }

    RepeatGenotype genotypeRepeatWithOneAllele(
        const CountTable& flanking_size_count, const CountTable& spanning_size_count,
        const std::vector<int32_t>& allele_size_candidates) const;

    RepeatGenotype genotypeRepeatWithTwoAlleles(
        const CountTable& flanking_size_count, const CountTable& spanning_size_count,
        const std::vector<int32_t>& allele_size_candidates) const;

private:
    int32_t repeatUnitLen_;
    int32_t max_repeat_size_in_units_;
    double prop_correct_molecules_;
};

}
