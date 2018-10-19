//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include <cstdint>
#include <memory>
#include <sstream>
#include <string>

#include "classification/alignment_classifier.h"
#include "graphalign/GraphAlignment.hh"

namespace reads
{

struct Read
{
    Read()
        : is_first_mate(false)
    {
    }

    Read(const std::string& new_read_id, const std::string& new_sequence)
        : read_id(new_read_id)
        , sequence(new_sequence)
        , is_first_mate(false)
    {
    }
    std::string read_id;
    std::string sequence;
    bool is_first_mate = false;

    std::string fragmentId() const
    {
        std::string fragment_id = read_id;
        fragment_id.pop_back();
        fragment_id.pop_back();
        return fragment_id;
    }

    const std::string& readId() const { return read_id; }
    bool isSecondMate() const { return !is_first_mate; }
    bool isSet() const { return !read_id.empty() && !sequence.empty(); }
};

struct LinearAlignmentStats
{
    int32_t chrom_id = -1;
    int32_t pos = -1;
    int32_t mapq = -1;
    int32_t mate_chrom_id = -1;
    int32_t mate_pos = -1;
    bool is_mapped = false;
    bool is_mate_mapped = false;
};

using ReadIdToLinearAlignmentStats = std::unordered_map<std::string, LinearAlignmentStats>;

bool operator==(const Read& read_a, const Read& read_b);
bool operator==(const LinearAlignmentStats& stats_a, const LinearAlignmentStats& core_info_b);

class RepeatAlignmentStats
{
public:
    RepeatAlignmentStats(
        const GraphAlignment& canonical_alignment, AlignmentType canonical_alignment_type,
        int32_t num_repeat_units_spanned)
        : canonical_alignment_(canonical_alignment)
        , canonical_alignment_type_(canonical_alignment_type)
        , num_repeat_units_spanned_(num_repeat_units_spanned)
    {
    }

    const GraphAlignment& canonicalAlignment() const { return canonical_alignment_; }
    AlignmentType canonicalAlignmentType() const { return canonical_alignment_type_; }
    int32_t numRepeatUnitsSpanned() const { return num_repeat_units_spanned_; }

private:
    GraphAlignment canonical_alignment_;
    AlignmentType canonical_alignment_type_;
    int32_t num_repeat_units_spanned_;
};

using ReadIdToRepeatAlignmentStats = std::unordered_map<std::string, RepeatAlignmentStats>;

std::ostream& operator<<(std::ostream& os, const Read& read);

} // namespace reads
