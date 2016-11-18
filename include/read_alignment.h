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

#ifndef INCLUDE_READ_ALIGNMENT_H_
#define INCLUDE_READ_ALIGNMENT_H_

#include <boost/lexical_cast.hpp>
using boost::lexical_cast;

#include <string>
#include <stdexcept>
#include <vector>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "include/genomic_region.h"

enum ReadStatus { kNoCheck, kFlankingRead };

struct Align {
  std::string name;
  int32_t len;
  std::string bases;
  std::string quals;
  int32_t chrom_id;
  int32_t pos;
  uint16_t mapq;
  uint32_t flag;
  int32_t mate_chrom_id;
  int32_t mate_pos;
  std::string region;
  ReadStatus status;

  bool IsMapped() const { return ((flag & 0x0004) == 0); }
  bool IsFirstMate() const { return ((flag & 0x0040) != 0); }
  bool IsMateMapped() const { return ((flag & 0x0008) == 0); }

  bool getMateRegion(Region& mateRegion,
                     const std::vector<std::string>& refVec) const {
    // Relies on mate (alignment) length being same as read (alignment) length.
    // May not hold for split alignments (or even gapped alignments?).
    if (flag & 0x08) {
      return false;  // mate unmapped
    }
    // align.mate_pos is 0-offset
    mateRegion = Region(DecodeChrom(mate_chrom_id, refVec), mate_pos + 1,
                        mate_pos + len);
    return true;
  }

  bool GetReadRegion(Region& readRegion,
                     const std::vector<std::string>& refVec) const {
    if (IsMapped()) {
      // myAlign.pos is 0-offset
      readRegion = Region(DecodeChrom(chrom_id, refVec), pos + 1, pos + len);
      return true;
    }

    return false;
  }

  std::string DecodeChrom(int32_t chromNum,
                          const std::vector<std::string>& refVec) const {
    if (chromNum == -1) {
      return "chr-1";
    }

    if (chromNum >= refVec.size()) {
      throw std::out_of_range(
          "[DecodeChrom ERROR] Input chromosme index: " +
          lexical_cast<std::string>(chromNum) + " but there are only " +
          lexical_cast<std::string>(refVec.size()) + " references");
    }

    return refVec[chromNum];
  }
};

bool GetAlignFromHtsAlign(bam1_t* hts_align_ptr, Align& align,
                          bool assumeUnaligned = false);

bool GetQualsFromHtsAlign(bam1_t* hts_align_ptr, std::string& quals);
bool GetBasesFromHtsAlign(bam1_t* hts_align_ptr, std::string& bases);

#endif  // INCLUDE_READ_ALIGNMENT_H_
