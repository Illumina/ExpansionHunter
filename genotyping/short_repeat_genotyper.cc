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

#include "genotyping/short_repeat_genotyper.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "common/common.h"

using std::string;
using std::map;
using std::cerr;
using std::endl;
using std::pair;
using std::vector;

Allele::Allele(int num_units_haplotype, int max_num_units_in_read,
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

double Allele::propMolecules(int num_units) const {
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

double Allele::propMoleculesShorterThan(int num_units_upper_bound) const {
  double total_prop = 0;
  for (int num_units = 0; num_units != num_units_upper_bound; ++num_units) {
    total_prop += propMolecules(num_units);
  }

  return total_prop;
}

double Allele::propMoleculesAtLeast(int num_units_lower_bound) const {

  return 1.0 - propMoleculesShorterThan(num_units_lower_bound);
}

double ShortRepeatGenotyper::CalcFlankingLoglik(int num_units_in_read) const {
  const double prob_start = hap_depth_ / read_len_;
  double loglik_sum = 0.0;
  for (const Allele &haplotype : alleles) {
    loglik_sum +=
        prob_start * haplotype.propMoleculesAtLeast(num_units_in_read);
  }

  const double gen_flanking_lik = loglik_sum / alleles.size();
  return log(gen_flanking_lik);
}

double ShortRepeatGenotyper::CalcSpanningLoglik(int num_units_in_read) const {
  const double prob_start = hap_depth_ / read_len_;
  double loglik_sum = 0.0;
  for (const Allele &haplotype : alleles) {
    loglik_sum += prob_start * haplotype.propMolecules(num_units_in_read);
  }
  const double gen_spanning_lik = loglik_sum / alleles.size();
  return log(gen_spanning_lik);
}

double
ShortRepeatGenotyper::CalcLogLik(const map<int, int> &flanking_size_counts,
                                 const map<int, int> &spanning_size_counts,
                                 std::vector<AlleleSupport> &support) const {
  double genotype_loglik = 0;
  support.resize(alleles.size());
  for (int i = 0; i != support.size(); ++i) {
    support[i] = {0, 0, 0};
  }

  for (const auto &kv : flanking_size_counts) {
    const int num_units = kv.first;
    int read_count = kv.second;

    int adjusted_read_count = read_count;
    if (num_units == max_num_units_in_read_) {
      adjusted_read_count = std::min<int>(read_count, 5);
    }

    genotype_loglik += adjusted_read_count * CalcFlankingLoglik(num_units);

    for (int i = 0; i != alleles.size(); ++i) {
      const int hap_num_units = alleles[i].num_units();
      if (num_units == max_num_units_in_read_) {
        if (hap_num_units == max_num_units_in_read_) {
          support[i].set_num_inrepeat(support[i].num_inrepeat() + read_count);
        }
      } else {
        if (num_units <= hap_num_units) {
          support[i].set_num_flanking(support[i].num_flanking() + read_count);
        }
      }
    }
  }

  for (const auto &kv : spanning_size_counts) {
    const int num_units = kv.first;
    const int read_count = kv.second;
    genotype_loglik += read_count * CalcSpanningLoglik(num_units);

    for (int i = 0; i != alleles.size(); ++i) {
      const int hap_num_units = alleles[i].num_units();
      if (num_units == hap_num_units) {
        support[i].set_num_spanning(support[i].num_spanning() + read_count);
      }
    }
  }

  return genotype_loglik;
}

void GenotypeShortRepeat(int max_num_units_in_read,
                         double prop_correct_molecules, double hap_depth,
                         int read_len,
                         const vector<RepeatAllele> &haplotype_candidates,
                         const map<int, int> &flanking_size_count,
                         const map<int, int> &spanning_size_count,
                         GenotypeType genotype_type, RepeatGenotype &genotype) {

  bool is_first_pass = true;
  double max_loglik = std::numeric_limits<double>::lowest();
  RepeatGenotype most_likely_genotype;

  if (genotype_type == GenotypeType::kDiploid) {
    for (const RepeatAllele &allele1 : haplotype_candidates) {
      for (const RepeatAllele &allele2 : haplotype_candidates) {
        if (allele1.size_ > allele2.size_) {
          continue;
        }
        ShortRepeatGenotyper genotype(max_num_units_in_read,
                                      prop_correct_molecules, hap_depth,
                                      read_len, allele1.size_, allele2.size_);
        vector<AlleleSupport> genotype_support;
        const double cur_loglik = genotype.CalcLogLik(
            flanking_size_count, spanning_size_count, genotype_support);

        if (max_loglik < cur_loglik || is_first_pass) {
          max_loglik = cur_loglik;
          is_first_pass = false;
          most_likely_genotype = {allele1, allele2};
          assert(most_likely_genotype.size() == genotype_support.size());
          most_likely_genotype[0].support_ = genotype_support[0];
          most_likely_genotype[1].support_ = genotype_support[1];
        }
      }
    }
  } else if (genotype_type == GenotypeType::kHaploid) {
    for (const RepeatAllele &allele : haplotype_candidates) {
      ShortRepeatGenotyper genotype(max_num_units_in_read,
                                    prop_correct_molecules, hap_depth, read_len,
                                    allele.size_);
      vector<AlleleSupport> genotype_support;
      const double cur_loglik = genotype.CalcLogLik(
          flanking_size_count, spanning_size_count, genotype_support);

      if (max_loglik < cur_loglik || is_first_pass) {
        max_loglik = cur_loglik;
        is_first_pass = false;
        most_likely_genotype = {allele};
        assert(most_likely_genotype.size() == genotype_support.size());
        most_likely_genotype[0].support_ = genotype_support[0];
      }
    }
  } else {
    throw std::logic_error("ERROR: Unknown GenotypeRepeat type");
  }
  genotype = most_likely_genotype;
};
