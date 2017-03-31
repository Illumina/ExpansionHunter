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

#ifndef INCLUDE_PARAMETERS_H_
#define INCLUDE_PARAMETERS_H_

#include <fstream>
#include <sstream>
#include <string>

#include "include/genomic_region.h"
#include "include/repeat_spec.h"
#include "genotyping/genotyping.h"

class Outputs {
 public:
  Outputs(const std::string vcf_path, const std::string json_path,
          const std::string log_path);
  std::ostream& vcf() { return vcf_; }
  std::ostream& json() { return json_; }
  std::ostream& log() { return log_; }

 private:
  std::ofstream vcf_;
  std::ofstream json_;
  std::ofstream log_;
};

class Parameters {
 public:
  const double kSmallestPossibleDepth = 0.00001;
  Parameters();
  bool Load(int numArgs, char* argPtrArr[]);
  const std::string& bam_path() const { return bam_path_; }
  const std::string& genome_path() const { return genome_path_; }
  const size_t region_extension_len() const { return region_extension_len_; }
  const float min_wp() const { return min_wp_; }
  const size_t min_baseq() const { return min_baseq_; }
  const size_t min_anchor_mapq() const { return min_anchor_mapq_; }
  const bool onlyUnaligned() const;
  const bool skip_unaligned() const { return skip_unaligned_; }
  const size_t read_len() const { return read_len_; };
  void set_read_len(size_t read_len) { read_len_ = read_len; }
  const double depth() const { return depth_; }
  void set_depth(double depth) { depth_ = depth; }
  const std::string& sample_name() const { return sample_name_; }
  const std::string& repeat_specs_path() const { return repeat_specs_path_; }
  const std::string& vcf_path() const { return vcf_path_; }
  const std::string& json_path() const { return json_path_; }
  const std::string& log_path() const { return log_path_; }
  bool depth_is_set() const { return depth_ >= kSmallestPossibleDepth; }

 private:
  std::string bam_path_;
  std::string genome_path_;
  // Specifies maximum distance from a target locus where
  // interesting reads may be.
  size_t region_extension_len_;
  float min_wp_;
  size_t min_baseq_;
  size_t min_anchor_mapq_;
  size_t read_len_;
  double depth_;
  Sex sex_;
  bool skip_unaligned_;
  std::string repeat_specs_path_;
  std::string sample_name_;
  std::string vcf_path_;
  std::string json_path_;
  std::string log_path_;
};

#endif  // INCLUDE_PARAMETERS_H_
