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

#include "reads/unaligned_reader.h"

#include <stdexcept>
#include <vector>

#include <boost/filesystem.hpp>

extern "C" {
#include "htslib/hts.h"
#include "htslib/sam.h"
}

using std::string;
using std::vector;

namespace reads {

struct UnalignedReader::Impl {
  Impl(const string &bam_path, const string &reference_path);
  Impl() = default;
  Impl(Impl const &) = delete;
  Impl(Impl &&) = delete;
  ~Impl();

  void PrepareBamForReading();
  void ExtractContigNamesFromHeader();
  void LoadIndex();

  void DestroyAlignmentPointer();
  void CloseRegion();
  void CloseFile();

  bool JumpToUnaligned();

  bool GetRead(Read &read);

  FileFormat file_format;
  string bam_path;
  string reference_path;

  vector<string> contig_names;
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
};

UnalignedReader::Impl::Impl(const string &new_bam_path,
                            const string &new_reference_path)
    : bam_path(new_bam_path), reference_path(new_reference_path) {
  PrepareBamForReading();
  ExtractContigNamesFromHeader();
  LoadIndex();
}

UnalignedReader::Impl::~Impl() {
  DestroyAlignmentPointer();
  CloseRegion();
  CloseFile();
}

void UnalignedReader::Impl::PrepareBamForReading() {
  hts_file_ptr = sam_open(bam_path.c_str(), "r");

  if (!hts_file_ptr) {
    throw std::runtime_error("Failed to read BAM/CRAM file " + bam_path);
  }

  switch (hts_file_ptr->format.format) {
    case bam:
      file_format = kBamFile;
      break;
    case cram:
      file_format = kCramFile;
    default:
      throw std::logic_error(
          "Could not identify file format of BAM/CRAM file " + bam_path);
  }

  if (file_format == kCramFile) {
    const string ref_index_path = reference_path + ".fai";
    if (!boost::filesystem::exists(ref_index_path)) {
      throw std::runtime_error("Reference index does not exist: " +
                               ref_index_path);
    }

    if (hts_set_fai_filename(hts_file_ptr, ref_index_path.c_str()) != 0) {
      throw std::runtime_error("Failed to set reference index");
    }
  }
}

void UnalignedReader::Impl::ExtractContigNamesFromHeader() {
  hts_bam_hdr_ptr = sam_hdr_read(hts_file_ptr);

  if (!hts_bam_hdr_ptr) {
    throw std::runtime_error("Failed to read header of: " + bam_path);
  }

  const int32_t num_contigs = hts_bam_hdr_ptr->n_targets;

  for (int32_t contig_index = 0; contig_index != num_contigs; ++contig_index) {
    const string contig_name(hts_bam_hdr_ptr->target_name[contig_index]);
    contig_names.push_back(contig_name);
  }
}

void BamReader::Impl::LoadIndex() {
  hts_idx_ptr = sam_index_load(hts_file_ptr, bam_path.c_str());

  if (!hts_idx_ptr) {
    throw std::runtime_error("Failed to read index of " + bam_path);
  }
}

void BamReader::Impl::DestroyAlignmentPointer() {
  if (hts_bam_align_ptr) {
    bam_destroy1(hts_bam_align_ptr);
    hts_bam_align_ptr = nullptr;
  }
}

void BamReader::Impl::CloseRegion() {
  if (hts_itr_ptr) {
    hts_itr_destroy(hts_itr_ptr);
    hts_itr_ptr = nullptr;
  }

  at_file_end = false;
}

void BamReader::Impl::CloseFile() {
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

  jump_to_unaligned = false;
  at_file_end = false;
}

}  // namespace reads