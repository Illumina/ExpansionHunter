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

#ifndef INCLUDE_REPEAT_SPEC_H_
#define INCLUDE_REPEAT_SPEC_H_

/*****************************************************************************/

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>

#include "include/genomic_region.h"

/*****************************************************************************/

class RepeatSpec {
 public:
  RepeatSpec() {}
  explicit RepeatSpec(const std::string& json_path);
  const char LeftFlankBase() const;
  bool is_common_unit() const { return is_common_unit_; }

  std::string repeat_id;
  Region target_region;
  std::string left_flank;
  std::string right_flank;
  std::string ref_seq;
  std::vector<std::string> units;
  std::vector<std::vector<std::string>> units_shifts;
  std::vector<Region> offtarget_regions;

 private:
  bool is_common_unit_;
};

bool LoadRepeatSpecs(const std::string& specs_path,
                     const std::string& genome_path, double min_wp,
                     std::map<std::string, RepeatSpec>* repeat_specs);

#endif  // INCLUDE_REPEAT_SPEC_H_
