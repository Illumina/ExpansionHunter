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

#include "genotyping/ShortRepeatGenotyper.hh"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "common/Common.hh"

using std::map;
using std::pair;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

QuantifierOfMoleculesGeneratedByAllele::QuantifierOfMoleculesGeneratedByAllele(
    int32_t allele_size_in_units, int32_t max_repeat_size_in_units, double prop_correct_molecules)
    : allele_size_in_units_(allele_size_in_units)
    , max_repeat_size_in_units_(max_repeat_size_in_units)
    , prop_correct_molecules_(prop_correct_molecules)
    , max_deviation_(5)
{
    norm_factor_ = 0;
    for (int num_units = 0; num_units <= max_repeat_size_in_units_; ++num_units)
    {
        const double deviation = std::abs(num_units - allele_size_in_units_);

        if (deviation < max_deviation_)
        {
            norm_factor_ += prop_correct_molecules_ * pow(1 - prop_correct_molecules_, deviation);
        }
        else
        {
            norm_factor_ += prop_correct_molecules_ * pow(1 - prop_correct_molecules_, max_deviation_);
        }
    }
}

double QuantifierOfMoleculesGeneratedByAllele::propMoleculesOfGivenSize(int size_in_units) const
{
    if (size_in_units < 0 || max_repeat_size_in_units_ < size_in_units)
    {
        throw std::logic_error(
            "size_in_units = " + to_string(size_in_units) + " is outside of allowed range (0,"
            + to_string(max_repeat_size_in_units_) + ")");
    }
    double deviation = std::abs(size_in_units - allele_size_in_units_);
    deviation = deviation < max_deviation_ ? deviation : max_deviation_;

    const double prop = prop_correct_molecules_ * pow(1 - prop_correct_molecules_, deviation);

    return prop / norm_factor_;
}

double QuantifierOfMoleculesGeneratedByAllele::propMoleculesShorterThan(int32_t size_upper_bound_in_units) const
{
    double total_prop = 0;
    for (int repeat_size_in_units = 0; repeat_size_in_units != size_upper_bound_in_units; ++repeat_size_in_units)
    {
        total_prop += propMoleculesOfGivenSize(repeat_size_in_units);
    }

    return total_prop;
}

double QuantifierOfMoleculesGeneratedByAllele::propMoleculesAtLeast(int32_t size_lower_bound_in_units) const
{
    return 1.0 - propMoleculesShorterThan(size_lower_bound_in_units);
}

double ShortRepeatGenotypeLikelihoodEstimator::CalcFlankingLoglik(int32_t num_units_in_read) const
{
    double loglik_sum = 0.0;
    for (const auto& allele_quantifier : allele_quantifiers_)
    {
        loglik_sum += allele_quantifier.propMoleculesAtLeast(num_units_in_read);
    }

    const double gen_flanking_lik = loglik_sum / allele_quantifiers_.size();
    return log(gen_flanking_lik);
}

double ShortRepeatGenotypeLikelihoodEstimator::CalcSpanningLoglik(int32_t num_units_in_read) const
{
    double loglik_sum = 0.0;
    for (const auto& allele_quantifier : allele_quantifiers_)
    {
        loglik_sum += allele_quantifier.propMoleculesOfGivenSize(num_units_in_read);
    }
    const double gen_spanning_lik = loglik_sum / allele_quantifiers_.size();
    return log(gen_spanning_lik);
}

double ShortRepeatGenotypeLikelihoodEstimator::CalcLogLik(
    const CountTable& counts_of_flanking_reads, const CountTable& counts_of_spanning_reads) const
{
    double genotype_loglik = 0;

    for (const auto& kv : counts_of_flanking_reads)
    {
        const int num_units = kv.first;
        int read_count = kv.second;

        int adjusted_read_count = read_count;
        if (num_units == max_repeat_size_in_units_)
        {
            adjusted_read_count = std::min<int32_t>(read_count, 5);
        }

        genotype_loglik += adjusted_read_count * CalcFlankingLoglik(num_units);
    }

    for (const auto& kv : counts_of_spanning_reads)
    {
        const int num_units = kv.first;
        const int read_count = kv.second;
        genotype_loglik += read_count * CalcSpanningLoglik(num_units);
    }

    return genotype_loglik;
}

RepeatGenotype ShortRepeatGenotyper::genotypeRepeatWithOneAllele(
    const CountTable& flanking_size_count, const CountTable& spanning_size_count,
    const vector<int32_t>& allele_size_candidates) const
{
    vector<int32_t> most_likely_genotype;
    double max_loglik = std::numeric_limits<double>::lowest();

    for (int32_t candidate_allele_size : allele_size_candidates)
    {
        ShortRepeatGenotypeLikelihoodEstimator likelihood_estimator(
            max_repeat_size_in_units_, prop_correct_molecules_, { candidate_allele_size });
        const double cur_loglik = likelihood_estimator.CalcLogLik(flanking_size_count, spanning_size_count);

        if (max_loglik < cur_loglik)
        {
            max_loglik = cur_loglik;
            most_likely_genotype = { candidate_allele_size };
        }
    }

    return RepeatGenotype(repeatUnitLen_, most_likely_genotype);
}

RepeatGenotype ShortRepeatGenotyper::genotypeRepeatWithTwoAlleles(
    const CountTable& flanking_size_count, const CountTable& spanning_size_count,
    const vector<int32_t>& allele_size_candidates) const
{
    vector<int32_t> most_likely_genotype;
    double max_loglik = std::numeric_limits<double>::lowest();

    for (int32_t candidate_allele_size_a : allele_size_candidates)
    {
        for (int32_t candidate_allele_size_b : allele_size_candidates)
        {
            if (candidate_allele_size_a > candidate_allele_size_b)
            {
                continue;
            }

            const vector<int32_t> candidate_genotype = { candidate_allele_size_a, candidate_allele_size_b };

            // ShortRepeatGenotyper genotyper(max_repeat_size_in_units_, prop_correct_molecules_);
            ShortRepeatGenotypeLikelihoodEstimator likelihood_estimator(
                max_repeat_size_in_units_, prop_correct_molecules_, candidate_genotype);
            const double cur_loglik = likelihood_estimator.CalcLogLik(flanking_size_count, spanning_size_count);

            if (max_loglik < cur_loglik)
            {
                max_loglik = cur_loglik;
                most_likely_genotype = candidate_genotype;
            }
        }
    }

    return RepeatGenotype(repeatUnitLen_, most_likely_genotype);
}

}
