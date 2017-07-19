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

#include <map>
#include <ostream>
#include <string>


enum class ReadType { kSpanning, kFlanking, kInrepeat, kOther };
const std::map<ReadType, std::string> kReadTypeToString = {
    {ReadType::kInrepeat, "INREPEAT"},
    {ReadType::kSpanning, "SPANNING"},
    {ReadType::kFlanking, "FLANKING"},
    {ReadType::kOther, "OTHER"}};

struct Read {
  std::string name;
  std::string bases;
  std::string quals;
};

class AlleleSupport {
public:
  AlleleSupport() : num_spanning_(0), num_flanking_(0), num_inrepeat_(0) {}
  AlleleSupport(int num_spanning, int num_flanking, int num_inrepeat)
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

  bool operator==(const AlleleSupport &rhs) const {
    return num_spanning_ == rhs.num_spanning_ &&
           num_flanking_ == rhs.num_flanking_ &&
           num_inrepeat_ == rhs.num_inrepeat_;
  }

private:
  int num_spanning_;
  int num_flanking_;
  int num_inrepeat_;
};

struct Interval {
  Interval() : lower_bound_(-1), upper_bound_(-1) {}
  int lower_bound_;
  int upper_bound_;
  bool operator==(const Interval &rhs) const {
    return lower_bound_ == rhs.lower_bound_ && upper_bound_ == rhs.upper_bound_;
  }
  std::string ToString() const {
    return std::to_string(lower_bound_) + "-" + std::to_string(upper_bound_);
  }
};

struct RepeatAllele {
  RepeatAllele(int size, int num_supporting_reads, ReadType type)
      : size_(size), num_supporting_reads_(num_supporting_reads), type_(type) {}
  RepeatAllele(int size, ReadType type, AlleleSupport support)
      : size_(size), type_(type), num_supporting_reads_(-1), support_(support) {}
  bool operator==(const RepeatAllele &rhs) const {
    return size_ == rhs.size_ && ci_ == rhs.ci_ && support_ == rhs.support_ &&
           type_ == rhs.type_ &&
           num_supporting_reads_ == rhs.num_supporting_reads_;
  }
  int size_;
  Interval ci_;
  AlleleSupport support_; // TODO: Rename to "consistent".
  int num_supporting_reads_;
  ReadType type_;
};

typedef std::vector<RepeatAllele> RepeatGenotype;