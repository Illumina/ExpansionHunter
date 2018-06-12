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
#include <string>

class Region {
 public:
  friend std::istream &operator>>(std::istream &istrm, Region &region);
  friend std::ostream &operator<<(std::ostream &ostrm, const Region &region);

  Region();
  Region(const std::string &chrom, int64_t start, int64_t end,
         const std::string &labelStr = std::string());
  Region(const std::string &rangeStr,
         const std::string &labelStr = std::string());

  bool is_set() const { return (chrom_ != "chr0"); }
  bool operator<(const Region &other_region) const;

  bool Overlaps(const Region &other_region) const;

  Region Extend(int extension_len) const;

  const std::string &chrom() const { return chrom_; }
  int64_t start() const { return start_; }
  int64_t end() const { return end_; }
  const std::string &label() const { return label_; }

  void set_start(int64_t start) { start_ = start; }
  void set_end(int64_t end) { end_ = end; }
  void set_label(const std::string &label) { label_ = label; }

  const std::string ToString() const;

 private:
  std::string chrom_;
  int64_t start_;
  int64_t end_;
  std::string label_;
};

std::istream &operator>>(std::istream &istrm, Region &region);
std::ostream &operator<<(std::ostream &ostrm, const Region &region);
