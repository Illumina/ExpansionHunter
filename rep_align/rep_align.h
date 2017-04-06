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

#ifndef REP_ALIGN_H
#define REP_ALIGN_H

#include <string>
#include <vector>

#include "common/common.h"
#include "common/parameters.h"
#include "common/repeat_spec.h"

struct RepeatAlign {
  enum class Type {
    kSpanning,
    kFlanking,
    kAnchored,
    kAlignedIrrPair,
    kUnalignedIrrPair,
    kUnalignedIrrSingleton
  };
  Read read;
  Read mate;
  size_t left_flank_len;
  size_t right_flank_len;
  Type type;
  size_t size;
};

size_t CountUnitsAtOffset(const std::vector<std::string> &units,
                          const std::string &bases, size_t offset);

size_t GetOffsetMostUnits(const std::vector<std::string> &units,
                          const std::string &bases, size_t *max_unit_count);

// Aligns read to left flank of the repeat; the alignment is considered valid
// if the raw wp score is 2 greater than the alignment of the same piece of the
// read against the repeat.
bool AlignLeftFlank(const std::vector<std::string> &units,
                    const std::string &left_flank, const std::string &bases,
                    const std::string &quals, size_t offset_most_units,
                    size_t min_baseq, double min_wp_score,
                    size_t *left_flank_len, double *left_flank_score);

// Counterpart of AlignLeftFlank for the right flank of the repeat.
bool AlignRightFlank(const std::vector<std::string> &units,
                     const std::string &right_flank, const std::string &bases,
                     const std::string &quals, size_t offset_most_units,
                     size_t min_baseq, double min_wp_score,
                     size_t *right_flank_len, double *right_flank_score);

// Tries to align the read to the repeat in forward orientation. Returns true if
// an alignment was found and populates rep_align with the alignment info.
bool IsSpanningOrFlankingRead(const Parameters &params,
                              const RepeatSpec &repeat_spec,
                              const std::string &bases,
                              const std::string &quals, RepeatAlign *rep_align);

// Tries to align the read in forward and reverse orientations.
bool IsSpanningOrFlankingReadRc(const Parameters &params,
                                const RepeatSpec &repeat_spec,
                                const std::string &bases,
                                const std::string &quals,
                                RepeatAlign *rep_align);

// Returns true if a flanking or a spanning read alignment was found; rep_align
// holds the actual repeat alignment info.
bool AlignRead(const Parameters &params, const RepeatSpec &repeat_spec,
               const std::string &bases, const std::string &quals,
               RepeatAlign *rep_align);

#endif // REP_ALIGN_H
