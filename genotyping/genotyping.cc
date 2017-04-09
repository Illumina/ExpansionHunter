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

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

using std::string;
using std::map;
using std::cerr;
using std::endl;
using std::pair;
using std::vector;
using std::array;

StrHaplotype::StrHaplotype(int num_units_haplotype, int max_num_units_in_read,
                           double prop_correct_molecules)
    : num_units_haplotype_(num_units_haplotype),
      max_num_units_in_read_(max_num_units_in_read),
      prop_correct_molecules_(prop_correct_molecules), max_deviation_(5) {
  const double p = prop_correct_molecules_;
  norm_factor_ = 0;
  for (int num_units = 0; num_units <= max_num_units_in_read_; ++num_units) {
    const double deviation = std::abs(num_units - num_units_haplotype_);

    if (deviation < max_deviation_) {
      norm_factor_ += p * pow(1 - p, deviation);
    } else {
      norm_factor_ += p * pow(1 - p, max_deviation_);
    }
  }
}

double StrHaplotype::propMolecules(int num_units) const {
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

double StrHaplotype::propMoleculesShorterThan(int num_units_upper_bound) const {
  double total_prop = 0;
  for (int num_units = 0; num_units != num_units_upper_bound; ++num_units) {
    total_prop += propMolecules(num_units);
  }

  return total_prop;
}

double StrHaplotype::propMoleculesAtLeast(int num_units_lower_bound) const {

  return 1.0 - propMoleculesShorterThan(num_units_lower_bound);
}

double StrGenotype::calcFlankingLoglik(int num_units_in_read) const {
  const double prob_start = hap_depth_ / read_len_;
  double loglik_sum = 0.0;
  for (const StrHaplotype &haplotype : haplotypes) {
    loglik_sum +=
        prob_start * haplotype.propMoleculesAtLeast(num_units_in_read);
  }

  const double gen_flanking_lik = loglik_sum / haplotypes.size();
  return log(gen_flanking_lik);
}

double StrGenotype::calcSpanningLoglik(int num_units_in_read) const {
  const double prob_start = hap_depth_ / read_len_;
  double loglik_sum = 0.0;
  for (const StrHaplotype &haplotype : haplotypes) {
    loglik_sum += prob_start * haplotype.propMolecules(num_units_in_read);
  }
  const double gen_spanning_lik = loglik_sum / haplotypes.size();
  return log(gen_spanning_lik);
}

double StrGenotype::calcLogLik(const map<int, int> &flanking_size_counts,
                               const map<int, int> &spanning_size_counts,
                               std::vector<std::array<int, 3>> &support) const {
  enum { kInrepeat = 0, kSpanning = 1, kFlanking = 2 };
  double genotype_loglik = 0;
  support.resize(haplotypes.size());
  for (int i = 0; i != support.size(); ++i) {
    support[i] = {0, 0, 0};
  }

  for (const auto &kv : flanking_size_counts) {
    const int num_units = kv.first;
    const int read_count = kv.second;
    genotype_loglik += read_count * calcFlankingLoglik(num_units);

    for (int i = 0; i != haplotypes.size(); ++i) {
      const int hap_num_units = haplotypes[i].num_units();
      if (num_units == max_num_units_in_read_) {
        if (hap_num_units == max_num_units_in_read_) {
          support[i][kInrepeat] += read_count;
        }
      } else {
        if (num_units <= hap_num_units) {
          support[i][kFlanking] += read_count;
        }
      }
    }
  }

  for (const auto &kv : spanning_size_counts) {
    const int num_units = kv.first;
    const int read_count = kv.second;
    genotype_loglik += read_count * calcSpanningLoglik(num_units);

    for (int i = 0; i != haplotypes.size(); ++i) {
      const int hap_num_units = haplotypes[i].num_units();
      if (num_units == hap_num_units) {
        support[i][kSpanning] += read_count;
      }
    }
  }

  return genotype_loglik;
}

vector<int> genotypeOneUnitStr(int max_num_units_in_read,
                               double prop_correct_molecules, double hap_depth,
                               int read_len,
                               const std::vector<int> &haplotype_candidates,
                               const map<int, int> &flanking_size_count,
                               const map<int, int> &spanning_size_count,
                               GenotypeType genotype_type) {
  bool is_first_pass = true;
  double max_loglik = std::numeric_limits<double>::lowest();
  vector<int> most_likely_genotype;
  vector<array<int, 3>> genotype_support;

  if (genotype_type == GenotypeType::kDiploid) {
    for (int num_units_hap1 : haplotype_candidates) {
      for (int num_units_hap2 : haplotype_candidates) {
        if (num_units_hap1 > num_units_hap2) {
          continue;
        }
        StrGenotype genotype(max_num_units_in_read, prop_correct_molecules,
                             hap_depth, read_len, num_units_hap1,
                             num_units_hap2);
        const double cur_loglik = genotype.calcLogLik(
            flanking_size_count, spanning_size_count, genotype_support);

        if (max_loglik < cur_loglik || is_first_pass) {
          max_loglik = cur_loglik;
          is_first_pass = false;
          most_likely_genotype = {num_units_hap1, num_units_hap2};
        }
      }
    }
  } else if (genotype_type == GenotypeType::kHaploid) {
    for (int num_units : haplotype_candidates) {
      StrGenotype genotype(max_num_units_in_read, prop_correct_molecules,
                           hap_depth, read_len, num_units);
      const double cur_loglik = genotype.calcLogLik(
          flanking_size_count, spanning_size_count, genotype_support);

      if (max_loglik < cur_loglik || is_first_pass) {
        max_loglik = cur_loglik;
        is_first_pass = false;
        most_likely_genotype = {num_units};
      }
    }
  } else {
    throw std::logic_error("ERROR: Unknown genotype type");
  }

  return most_likely_genotype;
};
