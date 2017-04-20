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

#include <iostream>
#include <stack>
#include <string>
#include <string>
#include <vector>

// Include BAM processing from samtools
#include "htslib/hts.h"
#include "htslib/sam.h"

#include "common/genomic_region.h"
#include "common/parameters.h"
#include "include/read_alignment.h"

class CramFile {
public:
  std::vector<int64_t> CountAlignedReads(const std::string &cram_path,
                                         int num_chroms);
  bool GetUnalignedRead(Align &align);

private:
  htsFile *file_ptr_;
  bam_hdr_t *header_ptr_;
  bam1_t *align_ptr_;
  bool found_unaligned_reads_;
};

class BamFile {
public:
  enum FileFormat { kBamFile, kCramFile, kUnknownFormat };

  BamFile();
  ~BamFile();

  void Init(const std::string &path, const std::string &reference);
  bool Close();
  bool SetRegionToRange(const Region &gRange);
  bool CloseRegion();
  bool JumpToUnaligned();
  double CalcMedianDepth(Parameters &parameters, size_t read_len);

  // If this returns false, there are no more reads in the set.
  // Before calling again, caller should open a new read set using. e.g.
  // SetRegionToRange(), JumpToUnaligned().
  bool GetRead(Align &align);

  // Can also be used for shadows since they have alignment pos of mate.
  bool GetAlignedMate(const Align &align, Align &mateBAlign);

  const std::vector<std::string> &ref_vec() const { return ref_vec_; }
  FileFormat format() const { return format_; }

private:
  bool GetUnalignedPrRead(Align &align);
  int GetNextGoodRead();

  enum { kSupplimentaryAlign = 0x800, kSecondaryAlign = 0x100 };

  CramFile cram_suppliment;
  // Contians chromosomes and their sizes.
  std::vector<std::string> ref_vec_;
  FileFormat format_;
  std::string path_;

  // Pointer to the input BAM/CRAM file itself.
  htsFile *hts_file_ptr_;
  // A pointer to BAM header.
  bam_hdr_t *hts_bam_hdr_ptr_;
  // A pointer to BAM index.
  hts_idx_t *hts_idx_ptr_;
  // Pointer to the target region.
  hts_itr_t *hts_itr_ptr_;
  // A pointer to an alignment in the BAM file.
  bam1_t *hts_bam_align_ptr_;

  bool jump_to_unaligned_;
  bool at_file_end_;
};
