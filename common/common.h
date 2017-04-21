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

#include <string>

struct Read {
  std::string name;
  std::string bases;
  std::string quals;
};

class HaplotypeSupport {
public:
  HaplotypeSupport() : num_spanning_(0), num_flanking_(0), num_inrepeat_(0) {}
  HaplotypeSupport(int num_spanning, int num_flanking, int num_inrepeat)
      : num_spanning_(num_spanning), num_flanking_(num_flanking),
        num_inrepeat_(num_inrepeat) {}

  int num_spanning() const { return num_spanning_; }
  int num_flanking() const { return num_flanking_; }
  int num_inrepeat() const { return num_inrepeat_; }

  void set_num_spanning(int num_spanning) { num_spanning_ = num_spanning; }
  void set_num_flanking(int num_flanking) { num_flanking_ = num_flanking; }
  void set_num_inrepeat(int num_inrepeat) { num_inrepeat_ = num_inrepeat; }

  std::string ToString() const {
    return std::to_string(num_spanning_) + "-" + std::to_string(num_flanking_) +
           "-" + std::to_string(num_inrepeat_);
  }

  bool operator==(const HaplotypeSupport &rhs) const {
    return num_spanning_ == rhs.num_spanning_ &&
           num_flanking_ == rhs.num_flanking_ &&
           num_inrepeat_ == rhs.num_inrepeat_;
  }

private:
  int num_spanning_;
  int num_flanking_;
  int num_inrepeat_;
};