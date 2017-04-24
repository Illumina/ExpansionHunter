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

double Genotype::CalcFlankingLoglik(int num_units_in_read) const {
  const double prob_start = hap_depth_ / read_len_;
  double loglik_sum = 0.0;
  for (const Allele &haplotype : alleles) {
    loglik_sum +=
        prob_start * haplotype.propMoleculesAtLeast(num_units_in_read);
  }

  const double gen_flanking_lik = loglik_sum / alleles.size();
  return log(gen_flanking_lik);
}

double Genotype::CalcSpanningLoglik(int num_units_in_read) const {
  const double prob_start = hap_depth_ / read_len_;
  double loglik_sum = 0.0;
  for (const Allele &haplotype : alleles) {
    loglik_sum += prob_start * haplotype.propMolecules(num_units_in_read);
  }
  const double gen_spanning_lik = loglik_sum / alleles.size();
  return log(gen_spanning_lik);
}

double Genotype::CalcLogLik(const map<int, int> &flanking_size_counts,
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

    if (num_units == max_num_units_in_read_) {
      read_count = std::min<int>(read_count, 5);
    }

    genotype_loglik += read_count * CalcFlankingLoglik(num_units);

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

vector<int> Genotype::ExtractAlleleSizes() const {
  vector<int> allele_sizes;
  for (const Allele &allele : alleles) {
    allele_sizes.push_back(allele.num_units());
  }
  return allele_sizes;
}

void GenotypeRepeat(int max_num_units_in_read, double prop_correct_molecules,
                    double hap_depth, int read_len,
                    const std::vector<int> &haplotype_candidates,
                    const map<int, int> &flanking_size_count,
                    const map<int, int> &spanning_size_count,
                    GenotypeType genotype_type, Genotype &genotype,
                    vector<string> &genotype_ci,
                    vector<AlleleSupport> &support) {

  bool is_first_pass = true;
  double max_loglik = std::numeric_limits<double>::lowest();
  Genotype most_likely_genotype;
  vector<AlleleSupport> most_likely_genotype_support;

  if (genotype_type == GenotypeType::kDiploid) {
    for (int num_units_hap1 : haplotype_candidates) {
      for (int num_units_hap2 : haplotype_candidates) {
        if (num_units_hap1 > num_units_hap2) {
          continue;
        }
        Genotype genotype(max_num_units_in_read, prop_correct_molecules,
                          hap_depth, read_len, num_units_hap1, num_units_hap2);
        vector<AlleleSupport> genotype_support;
        const double cur_loglik = genotype.CalcLogLik(
            flanking_size_count, spanning_size_count, genotype_support);

        if (max_loglik < cur_loglik || is_first_pass) {
          max_loglik = cur_loglik;
          is_first_pass = false;
          most_likely_genotype = genotype;
          most_likely_genotype_support = genotype_support;
        }
      }
    }
  } else if (genotype_type == GenotypeType::kHaploid) {
    for (int num_units : haplotype_candidates) {
      Genotype genotype(max_num_units_in_read, prop_correct_molecules,
                        hap_depth, read_len, num_units);
      vector<AlleleSupport> genotype_support;
      const double cur_loglik = genotype.CalcLogLik(
          flanking_size_count, spanning_size_count, genotype_support);

      if (max_loglik < cur_loglik || is_first_pass) {
        max_loglik = cur_loglik;
        is_first_pass = false;
        most_likely_genotype = genotype;
        most_likely_genotype_support = genotype_support;
      }
    }
  } else {
    throw std::logic_error("ERROR: Unknown GenotypeRepeat type");
  }
  genotype = most_likely_genotype;
  genotype_ci.clear();
  for (int size : genotype.ExtractAlleleSizes()) {
    genotype_ci.push_back(".");
  }
  support = most_likely_genotype_support;
};


/*
    int len_estimate = 0;
    int lower_bound = 0;
    int upper_bound = 0;

    // Haplotype depth should be twice as high because flanking reads
    // are coming from both flanks.
    EstimateRepeatLen(num_reads_from_unseen_allele, read_len, 2 * hap_depth,
                      len_estimate, lower_bound, upper_bound);

    // estimateRepeatLen adds read_len to size estimates so
    // we need to subtract it.
    len_estimate -= read_len;
    lower_bound -= read_len;
    upper_bound -= read_len;

    len_estimate = len_estimate / motif_len + longest_spanning + 1;
    lower_bound = lower_bound / motif_len + longest_spanning + 1;
    upper_bound = upper_bound / motif_len + longest_spanning + 1;

    // Repeat must be at least at long as the longest flanking read.
    lower_bound = std::max(lower_bound, longest_flanking);
    len_estimate = std::max(len_estimate, longest_flanking);
    upper_bound = std::max(upper_bound, longest_flanking);

    // Repeat estimated from flanking reads cannot be longer than the read
    // length.
    const int num_rep_in_read = read_len / motif_len;
    lower_bound = std::min(lower_bound, num_rep_in_read);
    len_estimate = std::min(len_estimate, num_rep_in_read);
    upper_bound = std::min(upper_bound, num_rep_in_read);

    // Make sure that size estimates have expected properties.
    if (!(lower_bound <= len_estimate && len_estimate <= upper_bound)) {
      cerr << "\t[Warning CoalesceFlankingReads: Unexpected size estimates. "
           << "Repeat size is " << len_estimate << " (LB=" << lower_bound
           << " UB=" << upper_bound << ")]" << endl;
    }
 */

/*
     cerr << "\t[Estimating repeat length from IRRs]" << endl;
    const double haplotype_depth = parameters.depth() / 2;
    EstimateRepeatLen(region_findings.num_irrs, parameters.read_len(),
                      haplotype_depth, repeat.size, repeat.size_ci_lower,
                      repeat.size_ci_upper);

 */