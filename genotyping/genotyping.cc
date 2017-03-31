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

#include "genotyping/genotyping.h"

#include <string>
#include <map>
#include <iostream>
#include <cctype>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

using std::string;
using std::map;
using std::cerr;
using std::endl;
using std::pair;
using std::vector;

StrHaplotype::StrHaplotype(int num_units_haplotype, int max_num_units_in_read, double prop_correct_molecules)
    : num_units_haplotype_(num_units_haplotype), max_num_units_in_read_(max_num_units_in_read),
      prop_correct_molecules_(prop_correct_molecules), max_deviation_(5)
{
  const double p = prop_correct_molecules_;
  norm_factor_ = 0;
  for (int num_units = 0; num_units <= max_num_units_in_read_; ++num_units) {
    const double deviation = std::abs(num_units - num_units_haplotype_);

    double prop = 0;
    if (deviation < max_deviation_) {
      norm_factor_ += p * pow(1 - p, deviation);
    } else {
      norm_factor_ += p * pow(1 - p, max_deviation_);
    }
  }
}

double StrHaplotype::propMolecules(int num_units)
{
  if (num_units < 0 || max_num_units_in_read_ < num_units) {
    cerr << "ERROR: num_units = " << num_units << endl;
  }
  const double p = prop_correct_molecules_;
  const double deviation = std::abs(num_units - num_units_haplotype_);

  double prop = 0;
  if (deviation < max_deviation_) {
    prop = p * pow(1 - p, deviation);
  } else {
    prop = p * pow(1 - p, max_deviation_);
  }

  return prop / norm_factor_;
}

double StrHaplotype::propMoleculesShorterThan(int num_units_upper_bound)
{
  double total_prop = 0;
  for (int num_units = 0; num_units != num_units_upper_bound; ++num_units) {
    total_prop += propMolecules(num_units);
  }

  return total_prop;
}

double StrHaplotype::propMoleculesAtLeast(int num_units_lower_bound) {

  return 1.0 - propMoleculesShorterThan(num_units_lower_bound);
}

double StrGenotype::calcFlankingLoglik(int num_units_in_read) {
  const double prob_start = hap_depth_ / read_len_;
  const double hap1_flanking_lik = prob_start * hap1_.propMoleculesAtLeast(num_units_in_read);
  const double hap2_flanking_lik = prob_start * hap2_.propMoleculesAtLeast(num_units_in_read);
  const double gen_flanking_lik = 0.5 * (hap1_flanking_lik + hap2_flanking_lik);
  return log(gen_flanking_lik);
}

double StrGenotype::calcSpanningLoglik(int num_units_in_read) {
  const double prob_start = hap_depth_ / read_len_;
  const double hap1_spanning_lik = prob_start * hap1_.propMolecules(num_units_in_read);
  const double hap2_spanning_lik = prob_start * hap2_.propMolecules(num_units_in_read);
  const double gen_spanning_lik = 0.5 * (hap1_spanning_lik + hap2_spanning_lik);
  return log(gen_spanning_lik);
}

double StrGenotype::calcLogLik(const map<int, int>& flanking_size_counts, const map<int, int>& spanning_size_counts)
{
  double genotype_loglik = 0;

  for (const auto& kv : flanking_size_counts) {
    const int num_units = kv.first;
    const int read_count = kv.second;
    genotype_loglik += read_count * calcFlankingLoglik(num_units);
  }

  for (const auto& kv : spanning_size_counts) {
    const int num_units = kv.first;
    const int read_count = kv.second;
    genotype_loglik += read_count * calcSpanningLoglik(num_units);
  }

  return genotype_loglik;
}

pair<int, int> genotypeOneUnitStr(int max_num_units_in_read, double prop_correct_molecules, double hap_depth,
                                  int read_len, const std::vector<int>& haplotype_candidates,
                                  const map<int, int>& flanking_size_count, const map<int, int>& spanning_size_count)
{
  bool is_first_pass = true;
  double max_loglik = std::numeric_limits<double>::lowest();
  pair<int, int> most_likely_genotype {0, 0};

  for (int num_units_hap1 : haplotype_candidates) {
    for (int num_units_hap2 : haplotype_candidates) {
      if (num_units_hap1 > num_units_hap2) {
        continue;
      }
      StrGenotype genotype(max_num_units_in_read, prop_correct_molecules, hap_depth, read_len, num_units_hap1,
                           num_units_hap2);
      const double cur_loglik = genotype.calcLogLik(flanking_size_count, spanning_size_count);

      cerr << num_units_hap1 << "/" << num_units_hap2 << "\t" << cur_loglik << endl;

      if (max_loglik < cur_loglik || is_first_pass) {
        max_loglik = cur_loglik;
        is_first_pass = false;
        most_likely_genotype.first = num_units_hap1;
        most_likely_genotype.second = num_units_hap2;
      }
    }
  }

  return most_likely_genotype;
};
