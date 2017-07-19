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

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/genomic_region.h"
#include "common/parameters.h"
#include "common/repeat_spec.h"
#include "include/bam_file.h"
#include "include/read_alignment.h"
#include "rep_align/rep_align.h"

// Represents fragment alignments as pairs of read alignments.
typedef std::array<Align, 2> AlignPair;

// Represents a collection fragments as a map from fragment names to fragment
// Aligns.
typedef std::unordered_map<std::string, AlignPair> AlignPairs;

enum WhatToCache { kCacheAll, kCacheIrr };

// Extract alignment pairs from a specified region.
void CacheReadsFromRegion(
    const Region &region, const WhatToCache whatToCache,
    const std::vector<std::vector<std::string>> &units_shifts,
    double min_wp_score, BamFile *bam_file, AlignPairs *align_pairs);

void CountAnchoredIrrs(
    const BamFile &bam_file, const Parameters &parameters,
    const Region &target_neighborhood, const RepeatSpec &repeat_spec,
    const std::unordered_set<std::string> &ontarget_frag_names,
    AlignPairs &align_pairs, int &num_anchored_irrs,
    const std::vector<std::vector<std::string>> &units_shifts,
    std::vector<RepeatAlign> *anchored_irrs);

void FillinMates(BamFile &bam_file, AlignPairs &align_pairs,
                 const std::vector<std::vector<std::string>> &units_shifts,
                 double min_wp_score,
                 const std::unordered_set<std::string> &ontarget_frag_names);

// Count the number of in-repeat reads stored in an AlignPairs object.
// A fragment is in-repeat if both of the reads fuzzy match to the repeat
// sequence.
bool CountUnalignedIrrs(
    BamFile &bam_file, const Parameters &parameters, int &numInRepeatReads,
    const std::vector<std::vector<std::string>> &units_shifts,
    std::vector<RepeatAlign> *irr_rep_aligns);

int CountAlignedIrr(const BamFile &bam_file, const Parameters &parameters,
                    const AlignPairs &align_pairs,
                    std::map<std::string, int> &num_irrs_per_offtarget_region,
                    const std::vector<std::vector<std::string>> &units_shifts,
                    std::vector<RepeatAlign> *irr_rep_aligns);
