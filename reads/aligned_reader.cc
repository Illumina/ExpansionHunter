//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>
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

#include "reads/aligned_reader.h"

#include <stdexcept>
#include <vector>

extern "C" {
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "common/genomic_region.h"
#include "reads/hts_helpers.h"

using std::string;
using std::vector;

namespace reads {

struct AlignedReader::Impl {
  Impl(const string &bam_path, const string &reference_path,
       const string &reference_index_path);
  Impl() = default;
  Impl(Impl const &) = delete;
  Impl(Impl &&) = delete;
  ~Impl();

  void PrepareBamForReading();
  void LoadHeader();
  void LoadIndex();

  void DestroyAlignmentPointer();
  void CloseRegion();
  void CloseFile();

  void SetRegion(const Region &region);
  bool GetRead(Read &read);
  int32_t ExtractNextPrimaryHtsAlignment();
  bool RecoverMate(const Read &read, Read &mate);

  string bam_path;
  string reference_path;
  string reference_index_path;

  // Pointer to the input BAM/CRAM file itself.
  htsFile *hts_file_ptr = nullptr;
  // A pointer to BAM header.
  bam_hdr_t *hts_bam_hdr_ptr = nullptr;
  // A pointer to BAM index.
  hts_idx_t *hts_idx_ptr = nullptr;
  // A pointer to the target region.
  hts_itr_t *hts_itr_ptr = nullptr;
  // A pointer to an alignment in the BAM file.
  bam1_t *hts_bam_align_ptr = nullptr;

  bool at_file_end = false;
};

AlignedReader::Impl::Impl(const string &new_bam_path,
                          const string &new_reference_path,
                          const string &reference_index_path)
    : bam_path(new_bam_path), reference_path(new_reference_path) {
  PrepareBamForReading();
  LoadHeader();
  LoadIndex();
}

AlignedReader::Impl::~Impl() {
  DestroyAlignmentPointer();
  CloseRegion();
  CloseFile();
}

void AlignedReader::Impl::PrepareBamForReading() {  // Open a BAM file for
                                                    // reading.
  hts_file_ptr = sam_open(bam_path.c_str(), "r");

  if (!hts_file_ptr) {
    throw std::runtime_error("Failed to read BAM/CRAM file " + bam_path);
  }

  if (hts_set_fai_filename(hts_file_ptr, reference_index_path.c_str()) != 0) {
    throw std::runtime_error("Failed to set reference index " +
                             reference_index_path);
  }
}

void AlignedReader::Impl::LoadHeader() {
  hts_bam_hdr_ptr = sam_hdr_read(hts_file_ptr);

  if (!hts_bam_hdr_ptr) {
    throw std::runtime_error("Failed to read header of: " + bam_path);
  }
}

void AlignedReader::Impl::LoadIndex() {
  hts_idx_ptr = sam_index_load(hts_file_ptr, bam_path.c_str());

  if (!hts_idx_ptr) {
    throw std::runtime_error("Failed to read index of " + bam_path);
  }
}

void AlignedReader::Impl::DestroyAlignmentPointer() {
  if (hts_bam_align_ptr) {
    bam_destroy1(hts_bam_align_ptr);
    hts_bam_align_ptr = nullptr;
  }
}

void AlignedReader::Impl::CloseRegion() {
  if (hts_itr_ptr) {
    hts_itr_destroy(hts_itr_ptr);
    hts_itr_ptr = nullptr;
  }

  at_file_end = false;
}

void AlignedReader::Impl::CloseFile() {
  if (hts_bam_hdr_ptr) {
    bam_hdr_destroy(hts_bam_hdr_ptr);
    hts_bam_hdr_ptr = nullptr;
  }

  if (hts_idx_ptr) {
    hts_idx_destroy(hts_idx_ptr);
    hts_idx_ptr = nullptr;
  }

  if (hts_file_ptr) {
    sam_close(hts_file_ptr);
    hts_file_ptr = nullptr;
  }

  at_file_end = false;
}

void AlignedReader::Impl::SetRegion(const Region &region) {
  // Close a currently open region, if any.
  CloseRegion();

  const string region_encoding = region.ToString();
  hts_itr_ptr =
      sam_itr_querys(hts_idx_ptr, hts_bam_hdr_ptr, region_encoding.c_str());

  if (hts_itr_ptr == 0) {
    throw std::runtime_error("Failed to extract reads from " + region_encoding);
  }

  at_file_end = false;
}

bool AlignedReader::Impl::GetRead(Read &read) {
  if (!hts_file_ptr) {
    throw std::logic_error("Cannot extract reads from closed BAM file");
  }

  if (!hts_itr_ptr) {
    throw std::logic_error("Read extraction requires target region to be set");
  }

  if (at_file_end) {
    return false;
  }

  if (!hts_bam_align_ptr) {
    hts_bam_align_ptr = bam_init1();
  }

  int32_t return_code = ExtractNextPrimaryHtsAlignment();

  if (return_code == -1) {
    at_file_end = true;
    return false;
  }

  if (return_code < -1) {
    throw std::runtime_error("Failed to extract read from BAM file");
  }

  htshelpers::DecodeAlignedRead(hts_bam_align_ptr, read);

  return true;
}

int32_t AlignedReader::Impl::ExtractNextPrimaryHtsAlignment() {
  while (true) {
    int32_t return_code =
        sam_itr_next(hts_file_ptr, hts_itr_ptr, hts_bam_align_ptr);
    if (return_code < 0) {  // Alignment extraction failed.
      return return_code;
    }
    const bool is_supplimentary =
        hts_bam_align_ptr->core.flag & htshelpers::kSupplementaryAlign;
    const bool is_secondary =
        hts_bam_align_ptr->core.flag & htshelpers::kSecondaryAlign;
    const bool is_primary_align = (!is_supplimentary) && (!is_secondary);
    if (is_primary_align) {
      return return_code;
    }
  }
}

bool AlignedReader::Impl::RecoverMate(const Read &read, Read &mate) {
  int32_t search_region_contig_id = 0;
  int32_t search_region_start = 0;
  int32_t search_region_end = 0;

  if (read.IsMateSamMapped()) {
    search_region_contig_id = read.SamMateChromId();
    search_region_start = read.SamMatePos();
  } else {
    search_region_contig_id = read.SamChromId();
    search_region_start = read.SamPos();
  }
  search_region_end = search_region_start + 1;

  hts_itr_t *hts_mate_region_iter_ptr;
  hts_mate_region_iter_ptr =
      sam_itr_queryi(hts_idx_ptr, search_region_contig_id, search_region_start,
                     search_region_end);

  if (!hts_mate_region_iter_ptr) {
    // cerr << "[Failed to parse " + mateRegion.ToString() << endl;
    return false;
  }

  // A pointer to an alignment in the BAM file.
  bam1_t *hts_possible_mate_align_ptr = nullptr;

  while (sam_itr_next(hts_file_ptr, hts_mate_region_iter_ptr,
                      hts_possible_mate_align_ptr) >= 0) {
    htshelpers::DecodeAlignedRead(hts_possible_mate_align_ptr, mate);

    if ((mate.FragmentId() == read.FragmentId()) &&
        (mate.IsFirstMate() != read.IsFirstMate())) {
      bam_destroy1(hts_possible_mate_align_ptr);
      hts_itr_destroy(hts_mate_region_iter_ptr);
      return true;
    }
  }
  bam_destroy1(hts_possible_mate_align_ptr);
  hts_itr_destroy(hts_mate_region_iter_ptr);
  return false;
}

}  // namespace reads