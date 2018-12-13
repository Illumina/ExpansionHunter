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

#include "common/common.h"
#include "common/parameters.h"
#include "common/repeat_spec.h"
#include "rep_align/rep_align.h"

namespace ehunter {

struct RepeatReadGroup {
  ReadType read_type;
  std::vector<RepeatAlign> rep_aligns;
  int size;
  int num_supporting_reads;
};

bool CompareReadGroupsBySize(const RepeatReadGroup &a1,
                             const RepeatReadGroup &a2);

void CoalesceFlankingReads(
    const RepeatSpec &repeat_spec, std::vector<RepeatReadGroup> &read_groups,
    std::vector<RepeatAlign> *flanking_repaligns, int motif_len,
    const std::vector<std::vector<std::string>> &units_shifts, int min_baseq,
    double min_wp_score);

void OutputRepeatAligns(const Parameters &parameters,
                        const RepeatSpec &repeat_spec,
                        const std::vector<RepeatReadGroup> &read_groups,
                        const std::vector<RepeatAlign> &flanking_repaligns,
                        std::ostream *out);

// Realign flanking reads to existing repeats.
void DistributeFlankingReads(const Parameters &parameters,
                             const RepeatSpec &repeat_spec,
                             std::vector<RepeatReadGroup> *read_groups,
                             std::vector<RepeatAlign> *flanking_repaligns);

}  // namespace ehunter
