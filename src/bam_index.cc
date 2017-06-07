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

#include "include/bam_index.h"

#include <boost/bind.hpp>

#include <iostream>
using std::cerr;
using std::endl;
#include <vector>
using std::vector;
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <cstdint>
#include <string>
using std::string;

#include "htslib/sam.h"

BamIndex::BamIndex(const string& bam_path) : bam_path_(bam_path) {
  // Set up hts_file_ptr_ (BAM/CRAM file pointer)
  hts_file_ptr_ = sam_open(bam_path_.c_str(), "r");

  if (!hts_file_ptr_) {
    cerr << "[ERROR]: Could not open '" << bam_path_ << "'" << endl;
    throw std::runtime_error("BamIndex: Could not open '" + bam_path_ + "'");
  }
}

BamIndex::~BamIndex() {
  sam_close(hts_file_ptr_);
  hts_file_ptr_ = 0;
}

bool BamIndex::GetChrReadCounts(vector<string>& chrom_names,
                                vector<int64_t>& chrom_lens,
                                vector<int64_t>& mapped_read_counts,
                                vector<int64_t>& unmapped_read_counts) const {
  chrom_names.clear();
  chrom_lens.clear();
  mapped_read_counts.clear();
  unmapped_read_counts.clear();

  bam_hdr_t* hts_bam_hdr_ptr = sam_hdr_read(hts_file_ptr_);

  if (hts_bam_hdr_ptr == 0) {
    cerr << "[ERROR]: GetChrReadCounts: "
         << "Failed to read header of BAM '" << bam_path_ << "'" << endl;
    return false;
  }

  hts_idx_t* hts_idx_ptr = sam_index_load(hts_file_ptr_, bam_path_.c_str());

  if (hts_idx_ptr == 0) {
    cerr << "[ERROR]: GetChrReadCounts : Failed to open index of BAM '"
         << bam_path_ << "'" << std::endl;
    return false;
  }

  const int chrom_count(hts_bam_hdr_ptr->n_targets);

  uint64_t mapped_count = 0;
  uint64_t unmapped_count = 0;

  for (int chrom_ind = 0; chrom_ind < chrom_count; ++chrom_ind) {
    chrom_names.push_back(string(hts_bam_hdr_ptr->target_name[chrom_ind]));
    chrom_lens.push_back(hts_bam_hdr_ptr->target_len[chrom_ind]);

    hts_idx_get_stat(hts_idx_ptr, chrom_ind, &mapped_count, &unmapped_count);

    mapped_read_counts.push_back((size_t)mapped_count);
    unmapped_read_counts.push_back((size_t)unmapped_count);
  }

  // Unaligned pairs
  chrom_names.push_back("*");
  chrom_lens.push_back(0);
  mapped_read_counts.push_back(0);
  unmapped_read_counts.push_back((size_t)hts_idx_get_n_no_coor(hts_idx_ptr));

  bam_hdr_destroy(hts_bam_hdr_ptr);
  hts_idx_destroy(hts_idx_ptr);

  return true;
}
