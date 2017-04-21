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

#pragma once

#include <boost/property_tree/ptree.hpp>

#include <array>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "common/parameters.h"
#include "common/repeat_spec.h"
#include "rep_align/rep_align.h"

struct Repeat {
  enum class SupportType { kInrepeat, kSpanning, kFlanking, kOther };
  std::map<SupportType, std::string> readtypeToStr = {
      {SupportType::kInrepeat, "INREPEAT"},
      {SupportType::kSpanning, "SPANNING"},
      {SupportType::kFlanking, "FLANKING"},
      {SupportType::kOther, "OTHER"}};
  std::vector<RepeatAlign> rep_aligns;
  int size;
  int size_ci_lower;
  int size_ci_upper;
  int num_supporting_reads;
  SupportType supported_by;
  void AsPtree(boost::property_tree::ptree &allele_node) const;
};

static bool CompareRepeatBySize(const Repeat &a1, const Repeat &a2) {
  return a1.size < a2.size;
}

// void AsPtree(const Parameters &parameters,
//             boost::property_tree::ptree &region_node,
//             std::vector<Repeat> alleles, const RepeatSpec &region_info,
//             const int num_irrs, const int num_unaligned_irrs,
//             const int num_anchored_irrs,
//             const std::vector<int> &off_target_irr_counts,
//             std::vector<int> &genotype,
//             const std::vector<std::array<int, 3>> &genotype_support);

void DumpVcf(const Parameters &parameters,
             const std::map<std::string, RepeatSpec> repeat_spec,
             const boost::property_tree::ptree &root_node, Outputs &outputs);

void CoalesceFlankingReads(
    const RepeatSpec &repeat_spec, std::vector<Repeat> &alleles,
    std::vector<RepeatAlign> *flanking_repaligns, const int read_len,
    const double hap_depth, int motif_len,
    const std::vector<std::vector<std::string>> &units_shifts, int min_baseq,
    double min_wp_score);

void OutputRepeatAligns(const Parameters &parameters,
                        const RepeatSpec &repeat_spec,
                        const std::vector<Repeat> &alleles,
                        const std::vector<RepeatAlign> &flanking_repaligns,
                        std::ostream *out);

// Realign flanking reads to existing repeats.
void DistributeFlankingReads(const Parameters &parameters,
                             const RepeatSpec &repeat_spec,
                             std::vector<Repeat> *alleles,
                             std::vector<RepeatAlign> *flanking_repaligns);
