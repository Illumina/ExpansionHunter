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

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "reads/read.h"

namespace reads {

class ReadPair {
 public:
  bool IsFirstMateSet() { return (bool)first_mate_ptr_; }
  bool IsSecondMateSet() { return (bool)second_mate_ptr_; }
  Read& FirstMate() { return *first_mate_ptr_; }
  Read& SecondMate() { return *second_mate_ptr_; }
  void Add(ReadPtr read_ptr);
  bool operator==(const ReadPair& other) const {
    return (first_mate_ptr_ == other.first_mate_ptr_ &&
            second_mate_ptr_ == other.second_mate_ptr_);
  }

  friend class ReadPairs;

 private:
  ReadPtr first_mate_ptr_;
  ReadPtr second_mate_ptr_;
};

/**
 * Read pair container class
 */
class ReadPairs {
 public:
  typedef std::unordered_map<std::string, ReadPair>::const_iterator
      const_iterator;
  const_iterator begin() const { return read_pairs_.begin(); }
  const_iterator end() const { return read_pairs_.end(); }

  ReadPairs() = default;
  void Clear();
  void Add(ReadPtr read_ptr);

  const ReadPair& operator[](const std::string& fragment_id) const;

  int32_t NumReads() const { return num_reads_; }

  void GetReads(std::vector<ReadPtr>& read_ptrs) const;

  bool operator==(const ReadPairs& other) const {
    return (read_pairs_ == other.read_pairs_ && num_reads_ == other.num_reads_);
  }

 private:
  std::unordered_map<std::string, ReadPair> read_pairs_;
  int32_t num_reads_ = 0;
};

}  // namespace reads