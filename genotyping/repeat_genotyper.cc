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

#include "genotyping/repeat_genotyper.h"

#include <cassert>
#include <iostream>

#include "genotyping/repeat_length.h"

using std::map;
using std::vector;
using std::string;
using std::cerr;
using std::endl;

void GenotypeRepeat(const Parameters &parameters, const RepeatSpec &repeat_spec,
                    int max_num_units_in_read, double prop_correct_molecules,
                    double hap_depth, int read_len,
                    const vector<RepeatAllele> &haplotype_candidates,
                    const map<int, int> &flanking_size_count,
                    const map<int, int> &spanning_size_count,
                    GenotypeType genotype_type, RepeatGenotype &genotype) {

  GenotypeShortRepeat(max_num_units_in_read, prop_correct_molecules, hap_depth,
                      read_len, haplotype_candidates, flanking_size_count,
                      spanning_size_count, genotype_type, genotype);

  const int unit_len = repeat_spec.units[0].length();
  const double haplotype_depth = parameters.depth() / 2;

  for (RepeatAllele& allele : genotype) {
    if (allele.type_ == ReadType::kSpanning) {
      allele.ci_.lower_bound_ = allele.size_;
      allele.ci_.upper_bound_ = allele.size_;
    }
  }

  if (genotype_type == GenotypeType::kDiploid) {
    RepeatAllele &short_allele = genotype.front();
    RepeatAllele &long_allele = genotype.back();
    assert(short_allele.size_ <= long_allele.size_);

    if (short_allele.type_ == ReadType::kInrepeat &&
        long_allele.type_ == ReadType::kInrepeat) {

      assert(flanking_size_count.at(max_num_units_in_read) ==
             long_allele.num_supporting_reads_);
      int num_irrs = long_allele.num_supporting_reads_;
      assert(num_irrs != 0);

      // Calculate CI for the short allele.
      int short_allele_size, short_allele_size_ci_lower,
          short_allele_size_ci_upper;

      EstimateRepeatLen(num_irrs / 2, parameters.read_len(), haplotype_depth,
                        short_allele_size, short_allele_size_ci_lower,
                        short_allele_size_ci_upper);

      short_allele_size /= unit_len;
      short_allele_size_ci_lower /= unit_len;
      short_allele_size_ci_upper /= unit_len;

      // Calculate CI for the long allele.
      int long_allele_size, long_allele_size_ci_lower,
          long_allele_size_ci_upper;

      EstimateRepeatLen(num_irrs, parameters.read_len(), haplotype_depth,
                        long_allele_size, long_allele_size_ci_lower,
                        long_allele_size_ci_upper);

      long_allele_size /= unit_len;
      long_allele_size_ci_lower /= unit_len;
      long_allele_size_ci_upper /= unit_len;

      short_allele.size_ = short_allele_size;
      short_allele.ci_.lower_bound_ = short_allele_size_ci_lower;
      short_allele.ci_.upper_bound_ = short_allele_size_ci_upper;

      long_allele.size_ = long_allele_size;
      long_allele.ci_.lower_bound_ = long_allele_size_ci_lower;
      long_allele.ci_.upper_bound_ = long_allele_size_ci_upper;
    } else if (long_allele.type_ == ReadType::kInrepeat) {
      assert(short_allele.type_ == ReadType::kSpanning);

      int num_irrs = flanking_size_count.at(max_num_units_in_read);
      assert(num_irrs != 0);
      int long_allele_size, long_allele_size_ci_lower,
          long_allele_size_ci_upper;
      EstimateRepeatLen(num_irrs, parameters.read_len(), haplotype_depth,
                        long_allele_size, long_allele_size_ci_lower,
                        long_allele_size_ci_upper);

      long_allele_size /= unit_len;
      long_allele_size_ci_lower /= unit_len;
      long_allele_size_ci_upper /= unit_len;

      long_allele.size_ = long_allele_size;
      long_allele.ci_.lower_bound_ = long_allele_size_ci_lower;
      long_allele.ci_.upper_bound_ = long_allele_size_ci_upper;
    } else if (long_allele.type_ == ReadType::kFlanking) {
      assert(short_allele.type_ == ReadType::kSpanning ||
             short_allele.type_ == ReadType::kFlanking);
      const int longest_flanking = long_allele.size_;
      int longest_spanning = 0;
      if (short_allele.type_ == ReadType::kSpanning) {
        longest_spanning = short_allele.size_;
      }
      int len_estimate = 0;
      int lower_bound = 0;
      int upper_bound = 0;
      // Haplotype depth should be twice as high because flanking reads
      // are coming from both flanks.
      EstimateRepeatLen(long_allele.num_supporting_reads_, read_len,
                        2 * hap_depth, len_estimate, lower_bound, upper_bound);

      // estimateRepeatLen adds read_len to size estimates so we need to
      // subtract it.
      len_estimate -= read_len;
      lower_bound -= read_len;
      upper_bound -= read_len;

      len_estimate = len_estimate / unit_len + longest_spanning + 1;
      lower_bound = lower_bound / unit_len + longest_spanning + 1;
      upper_bound = upper_bound / unit_len + longest_spanning + 1;

      // Repeat must be at least at long as the longest flanking read.
      lower_bound = std::max(lower_bound, longest_flanking);
      len_estimate = std::max(len_estimate, longest_flanking);
      upper_bound = std::max(upper_bound, longest_flanking);

      // Repeat estimated from flanking reads cannot be longer than the read
      // length.
      const int num_rep_in_read = read_len / unit_len;
      lower_bound = std::min(lower_bound, num_rep_in_read);
      len_estimate = std::min(len_estimate, num_rep_in_read);
      upper_bound = std::min(upper_bound, num_rep_in_read);

      long_allele.size_ = len_estimate;
      long_allele.ci_.lower_bound_ = lower_bound;
      long_allele.ci_.upper_bound_ = upper_bound;

      if (short_allele.type_ == ReadType::kFlanking) {
        short_allele.size_ = long_allele.size_;
        short_allele.ci_.lower_bound_ = long_allele.ci_.lower_bound_;
        short_allele.ci_.upper_bound_ = long_allele.ci_.upper_bound_;
      }

    } else {
      // Both alleles must be spanning.
      assert(short_allele.type_ == ReadType::kSpanning &&
             long_allele.type_ == ReadType::kSpanning);
    }
  } else {
    assert(genotype_type == GenotypeType::kHaploid);
    RepeatAllele &allele = genotype.front();
    if (allele.type_ == ReadType::kInrepeat) {
      int num_irrs = allele.num_supporting_reads_;
      assert(num_irrs != 0);
      int allele_size, allele_size_ci_lower, allele_size_ci_upper;
      EstimateRepeatLen(num_irrs, parameters.read_len(), haplotype_depth,
                        allele_size, allele_size_ci_lower,
                        allele_size_ci_upper);

      allele_size /= unit_len;
      allele_size_ci_lower /= unit_len;
      allele_size_ci_upper /= unit_len;

      allele.size_ = allele_size;
      allele.ci_.lower_bound_ = allele_size_ci_lower;
      allele.ci_.upper_bound_ = allele_size_ci_upper;
    } else if (allele.type_ == ReadType::kFlanking) {
      const int longest_flanking = allele.size_;
      int longest_spanning = 0;
      int len_estimate = 0;
      int lower_bound = 0;
      int upper_bound = 0;
      // Haplotype depth should be twice as high because flanking reads
      // are coming from both flanks.
      EstimateRepeatLen(allele.num_supporting_reads_, read_len,
                        2 * hap_depth, len_estimate, lower_bound, upper_bound);

      // estimateRepeatLen adds read_len to size estimates so we need to
      // subtract it.
      len_estimate -= read_len;
      lower_bound -= read_len;
      upper_bound -= read_len;

      len_estimate = len_estimate / unit_len + longest_spanning + 1;
      lower_bound = lower_bound / unit_len + longest_spanning + 1;
      upper_bound = upper_bound / unit_len + longest_spanning + 1;

      // Repeat must be at least at long as the longest flanking read.
      lower_bound = std::max(lower_bound, longest_flanking);
      len_estimate = std::max(len_estimate, longest_flanking);
      upper_bound = std::max(upper_bound, longest_flanking);

      // Repeat estimated from flanking reads cannot be longer than the read
      // length.
      const int num_rep_in_read = read_len / unit_len;
      lower_bound = std::min(lower_bound, num_rep_in_read);
      len_estimate = std::min(len_estimate, num_rep_in_read);
      upper_bound = std::min(upper_bound, num_rep_in_read);

      allele.size_ = len_estimate;
      allele.ci_.lower_bound_ = lower_bound;
      allele.ci_.upper_bound_ = upper_bound;
    }
  }
}
