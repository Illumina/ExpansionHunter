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

#include "include/bam_file.h"

#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/sam.h"

#include "common/parameters.h"
#include "common/ref_genome.h"
#include "common/timestamp.h"
#include "include/bam_index.h"

typedef boost::tokenizer<boost::char_separator<char>> Tokenizer;
using boost::lexical_cast;

using std::string;
using std::cerr;
using std::endl;
using std::vector;
using std::sort;

BamFile::BamFile()
    : hts_file_ptr_(0), hts_bam_hdr_ptr_(0), hts_idx_ptr_(0), hts_itr_ptr_(0),
      hts_bam_align_ptr_(0), jump_to_unaligned_(false), at_file_end_(false),
      format_(kUnknownFormat) {}

BamFile::~BamFile() {
  if (hts_bam_align_ptr_) {
    bam_destroy1(hts_bam_align_ptr_);
    hts_bam_align_ptr_ = 0;
  }

  Close();
}

void BamFile::Init(const string &path, const string &reference) {
  path_ = path;
  // Open a BAM file for reading.
  hts_file_ptr_ = sam_open(path.c_str(), "r");

  if (!hts_file_ptr_) {
    throw std::runtime_error("BamFile::Init: Failed to read BAM file '" + path +
                             "'");
  }

  if (hts_file_ptr_->format.format == bam) {
    format_ = kBamFile;
  }

  if (hts_file_ptr_->format.format == cram) {
    format_ = kCramFile;

    const string reference_index = reference + ".fai";
    if (!boost::filesystem::exists(reference_index)) {
      throw std::runtime_error("Reference index does not exist: " +
                               reference_index);
    }

    if (hts_set_fai_filename(hts_file_ptr_, reference_index.c_str()) != 0) {
      throw std::runtime_error("Failed to set reference index");
    }
  }

  string input_format = "Unknown";
  if (format_ == kBamFile) {
    input_format = "BAM";
  } else if (format_ == kCramFile) {
    input_format = "CRAM";
  }

  cerr << TimeStamp() << ",[Input format: " << input_format << "]" << endl;

  // Read hdr and set up ref_vec_
  hts_bam_hdr_ptr_ = sam_hdr_read(hts_file_ptr_);

  if (!hts_bam_hdr_ptr_) {
    throw std::runtime_error("BamFile::Init: Failed to read BAM header: '" +
                             path + "'");
  }

  const int chrom_count = hts_bam_hdr_ptr_->n_targets;

  for (int chrom_ind = 0; chrom_ind < chrom_count; ++chrom_ind) {
    ref_vec_.push_back(string(hts_bam_hdr_ptr_->target_name[chrom_ind]));
  }

  // Load the index
  hts_idx_ptr_ = sam_index_load(hts_file_ptr_, path.c_str());

  if (!hts_idx_ptr_) {
    throw std::runtime_error("BamFile::Init: Failed to read BAM index '" +
                             path + "'");
  }
}

bool BamFile::Close() {
  if (hts_file_ptr_) {
    if (hts_bam_hdr_ptr_) {
      bam_hdr_destroy(hts_bam_hdr_ptr_);
      hts_bam_hdr_ptr_ = 0;
    }

    CloseRegion();

    if (hts_idx_ptr_) {
      hts_idx_destroy(hts_idx_ptr_);
      hts_idx_ptr_ = 0;
    }

    sam_close(hts_file_ptr_);
    hts_file_ptr_ = 0;
  }

  jump_to_unaligned_ = false;
  at_file_end_ = false;

  return true;
}

// Set bam file to a specific range from which the reads will be extracted.
bool BamFile::SetRegionToRange(const Region &gRange) {
  // If we were in unaligned pairs reset first.
  if (jump_to_unaligned_) {
    jump_to_unaligned_ = false;
  }

  // If we were in another region, close that.
  if (hts_itr_ptr_ != 0) {
    CloseRegion();
  }

  hts_itr_ptr_ =
      sam_itr_querys(hts_idx_ptr_, hts_bam_hdr_ptr_, gRange.ToString().c_str());

  if (hts_itr_ptr_ == 0) {
    throw std::runtime_error("Failed to set target region: '" +
                             gRange.ToString() + "'");
    return false;
  }

  at_file_end_ = false;

  return true;
}

bool BamFile::CloseRegion() {
  if (hts_itr_ptr_) {
    hts_itr_destroy(hts_itr_ptr_);
    hts_itr_ptr_ = 0;
  }

  at_file_end_ = false;

  return true;
}

// Hack to get to the unaligned pairs (at the end of the BAM) quickly.
// Find the offset of the last block containing aligned reads, then skip
// reads until reach unaligned ones.
// Based on Isis UnalignedReadExtractor::JumpToUnalignedReads.

bool BamFile::JumpToUnaligned() {
  // If we were in another region, close that.
  if (hts_itr_ptr_ != 0) {
    CloseRegion();
    hts_itr_ptr_ = 0;
  }

  if (hts_file_ptr_->format.format == bam) {
    const bool hasUnalignedPairs(hts_idx_get_n_no_coor(hts_idx_ptr_) > 0);

    if (!hasUnalignedPairs) {
      return false;
    }

    hts_itr_ptr_ = sam_itr_querys(hts_idx_ptr_, hts_bam_hdr_ptr_, "*");

    if (hts_itr_ptr_ == 0) {
      throw std::runtime_error("Failed to extract an unaligned read");
    }

    jump_to_unaligned_ = true;
  } else if (hts_file_ptr_->format.format == cram) {
    jump_to_unaligned_ = true;
  } else {
    throw std::logic_error("Unknown format");
  }

  return true;
}

bool BamFile::GetRead(Align &align) {
  if (jump_to_unaligned_) {
    if (hts_file_ptr_->format.format == bam) {
      return GetUnalignedPrRead(align);
    } else if (hts_file_ptr_->format.format == cram) {
      return cram_suppliment.GetUnalignedRead(align);
    } else {
      throw std::logic_error("Unknown format");
    }
  }

  if (!hts_file_ptr_) {
    throw std::logic_error("BamFile::GetRead BAM file is not open");
  }

  if (at_file_end_) {
    return false;
  }

  if (!hts_bam_align_ptr_) {
    hts_bam_align_ptr_ = bam_init1();
  }

  int readRet = 0;
  assert(hts_itr_ptr_);
  readRet = GetNextGoodRead();

  if (readRet == -1) {
    at_file_end_ = true; // EOF
    return false;
  }

  if (readRet < -1) {
    throw std::runtime_error("Failed to extract read from BAM file");
  }

  if (!GetAlignFromHtsAlign(hts_bam_align_ptr_, align)) {
    throw std::runtime_error("Failed to process read from BAM file");
  }

  return true;
}

// Try to get an aligned mate.
bool BamFile::GetAlignedMate(const Align &align, Align &mate_align) {
  Region mateRegion;
  int32_t tid = 0;
  int32_t beg = 0;
  int32_t end = 0;

  if (align.IsMateMapped()) {
    tid = align.mate_chrom_id;
    beg = align.mate_pos;
    end = align.mate_pos + 1;
  } else {
    tid = align.chrom_id;
    beg = align.pos;
    end = align.pos + 1;
  }

  hts_itr_t *iter;
  iter = sam_itr_queryi(hts_idx_ptr_, tid, beg, end);
  if (!iter) {
    cerr << "[Failed to parse " + mateRegion.ToString() << endl;
    return false;
  }
  while (sam_itr_next(hts_file_ptr_, iter, hts_bam_align_ptr_) >= 0) {
    if (!GetAlignFromHtsAlign(hts_bam_align_ptr_, mate_align)) {
      hts_itr_destroy(iter);
      throw std::runtime_error("Failed to process read from BAM file");
    }
    if ((mate_align.name == align.name) &&
        (mate_align.IsFirstMate() != align.IsFirstMate())) {
      hts_itr_destroy(iter);
      return true;
    }
  }
  hts_itr_destroy(iter);
  return false;
}

bool BamFile::GetUnalignedPrRead(Align &align) {
  if (at_file_end_) {
    return false;
  }
  if (!hts_itr_ptr_) {
    throw std::runtime_error("GetUnalignedPrRead but hts_itr_ptr_ 0");
  }
  if (!hts_bam_align_ptr_) {
    hts_bam_align_ptr_ = bam_init1();
  }

  bool found_unaligned = false;

  // Skip any aligned reads.
  while (sam_itr_next(hts_file_ptr_, hts_itr_ptr_, hts_bam_align_ptr_) > 0) {
    // unaligned pair -> 0x4 unmapped + 0x8 mate unmapped -> 0xC
    if ((hts_bam_align_ptr_->core.flag & 0xC) == 0xC) {
      found_unaligned = true;
      break;
    }
  }

  if (!found_unaligned) {
    at_file_end_ = true;
    return false;
  }

  // Have retrieved an unaligned read. Copy the bits needed to align.
  if (!GetAlignFromHtsAlign(hts_bam_align_ptr_, align,
                            true)) { // assumeUnaligned=T
    throw std::runtime_error("Failed to process read from BAM file.");
  }

  return true;
}

int BamFile::GetNextGoodRead() {
  bool is_primary_align = false;
  int return_value = 0;

  while (!is_primary_align) {
    return_value =
        sam_itr_next(hts_file_ptr_, hts_itr_ptr_, hts_bam_align_ptr_);
    if (return_value < 0) {
      // low-level reading failed so report the return code.
      return return_value;
    }
    const bool is_supplimentary =
        hts_bam_align_ptr_->core.flag & kSupplimentaryAlign;
    const bool is_secondary = hts_bam_align_ptr_->core.flag & kSecondaryAlign;
    is_primary_align = (!is_supplimentary) && (!is_secondary);
  }
  return return_value;
}

const size_t CountValidBases(const string &bases) {
  const size_t n_count = std::count(bases.begin(), bases.end(), 'N');
  const size_t valid_base_count = bases.length() - n_count;

  return valid_base_count;
}

static bool IsAutosome(const string &chrom_name) {
  static std::unordered_set<string> autosome_names = {
      "chr1",  "chr2",  "chr3",  "chr4",  "chr5",  "chr6",  "chr7",  "chr8",
      "chr9",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
      "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "1",     "2",
      "3",     "4",     "5",     "6",     "7",     "8",     "9",     "10",
      "11",    "12",    "13",    "14",    "15",    "16",    "17",    "18",
      "19",    "20",    "21",    "22"};

  auto search = autosome_names.find(chrom_name);
  if (search != autosome_names.end()) {
    return true;
  }
  return false;
}

double BamFile::CalcMedianDepth(Parameters &parameters, size_t read_len) {
  if (read_len == 0) {
    throw std::logic_error("Read length must be non-zero: " +
                           lexical_cast<string>(read_len));
  }

  RefGenome ref_genome(parameters.genome_path());
  BamIndex bam_index(parameters.bam_path());

  vector<int64_t> mapped_read_counts;
  vector<int64_t> unmapped_read_counts;

  vector<string> chrom_names;
  vector<int64_t> chrom_lens;
  if (!bam_index.GetChrReadCounts(chrom_names, chrom_lens, mapped_read_counts,
                                  unmapped_read_counts)) {
    throw std::runtime_error("Failed to get chrom read depths from index of " +
                             parameters.bam_path());
  }

  const int chrom_count = chrom_names.size();

  if (format_ == kCramFile) {
    mapped_read_counts =
        cram_suppliment.CountAlignedReads(parameters.bam_path(), chrom_count);
    for (int i = 0; i < chrom_names.size(); ++i) {
      cerr << chrom_names[i] << " " << chrom_lens[i] << " "
           << mapped_read_counts[i] << endl;
    }
  }

  typedef std::pair<int, double> ChromIndDepth;
  typedef vector<ChromIndDepth> ChromIndDepths;
  ChromIndDepths chrom_ind_depths;

  string chrom_bases;

  for (int chrom_ind = 0; chrom_ind < chrom_count; ++chrom_ind) {
    const string &chrom_name = chrom_names[chrom_ind];

    if (!IsAutosome(chrom_name)) {
      continue;
    }

    cerr << TimeStamp() << ",[Using " << chrom_name << " to calculate depth]"
         << endl;

    ref_genome.ExtractSeq(chrom_name, &chrom_bases);

    const double read_depth =
        static_cast<double>(mapped_read_counts[chrom_ind]) * read_len /
        CountValidBases(chrom_bases);

    chrom_ind_depths.push_back(ChromIndDepth(chrom_ind, read_depth));
  }

  if (chrom_ind_depths.empty()) {
    throw std::runtime_error("Error: No contigs named chr1-chr22 or 1-22 "
                             "found; consider setting the depth manually");
  }

  sort(chrom_ind_depths.begin(), chrom_ind_depths.end(),
       boost::bind(&std::pair<int, double>::second, _1) >
           boost::bind(&std::pair<int, double>::second, _2));

  const size_t autosome_count = chrom_ind_depths.size();
  const bool autosome_count_is_odd = (autosome_count % 2) == 1;
  const size_t half_autosome_count = autosome_count / 2;

  const double median_autosome_depth =
      (autosome_count_is_odd
           ? (chrom_ind_depths[half_autosome_count].second)
           : ((chrom_ind_depths[half_autosome_count - 1].second +
               chrom_ind_depths[half_autosome_count].second) /
              2.0));
  return median_autosome_depth;
}

vector<int64_t> CramFile::CountAlignedReads(const string &cram_path,
                                            int num_chroms) {
  vector<int64_t> read_counts(num_chroms, 0);
  file_ptr_ = sam_open(cram_path.c_str(), "r");

  if (!file_ptr_) {
    throw std::runtime_error("Failed to read the input file");
  }

  if (file_ptr_->format.format != cram) {
    throw std::runtime_error(cram_path + " is not a CRAM file");
  }

  header_ptr_ = sam_hdr_read(file_ptr_);
  if (!header_ptr_) {
    throw std::runtime_error("Could not read header of " + cram_path);
  }

  align_ptr_ = bam_init1();

  found_unaligned_reads_ = false;
  int ret;
  while ((ret = sam_read1(file_ptr_, header_ptr_, align_ptr_)) >= 0) {
    if (align_ptr_->core.tid == -1) { // Reached unaligned reads.
      found_unaligned_reads_ = true;
      cerr << "[Found unaligend reads]" << endl;
      break;
    }
    // Counting mapped reads only.
    if ((align_ptr_->core.flag & 0x0004) == 0) {
      ++read_counts[align_ptr_->core.tid];
    }
  }

  return read_counts;
}

bool CramFile::GetUnalignedRead(Align &align) {
  int ret = sam_read1(file_ptr_, header_ptr_, align_ptr_);
  if (ret < 0) {
    return false;
  }
  if (!GetAlignFromHtsAlign(align_ptr_, align, true)) {
    throw std::runtime_error("Failed to process read from BAM file.");
  }
  return true;
}
