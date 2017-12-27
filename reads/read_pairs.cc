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

#include "reads/read_pairs.h"

#include <stdexcept>

using std::string;
using std::vector;

namespace reads {

void ReadPair::Add(ReadPtr read_ptr) {
  if (read_ptr->IsFirstMate()) {
    first_mate_ptr_ = read_ptr;
  } else {
    second_mate_ptr_ = read_ptr;
  }
}

void ReadPairs::Add(ReadPtr read_ptr) {
  ReadPair& read_pair = read_pairs_[read_ptr->FragmentId()];
  const int32_t num_mates_original = (int32_t)read_pair.IsFirstMateSet() +
                                     (int32_t)read_pair.IsSecondMateSet();
  read_pair.Add(read_ptr);
  const int32_t num_mates_after_add = (int32_t)read_pair.IsFirstMateSet() +
                                      (int32_t)read_pair.IsSecondMateSet();

  num_reads_ += num_mates_original - num_mates_after_add;
}

const ReadPair& ReadPairs::operator[](const string& fragment_id) const {
  if (read_pairs_.find(fragment_id) == read_pairs_.end()) {
    throw std::logic_error("Fragment " + fragment_id + " does not exist");
  }
  return read_pairs_.at(fragment_id);
}

void ReadPairs::GetReads(vector<ReadPtr>& reads) const {
  for (const auto& kv : read_pairs_) {
    ReadPair read_pair = kv.second;
    if (read_pair.IsFirstMateSet()) {
      reads.push_back(read_pair.first_mate_ptr_);
    }
    if (read_pair.IsSecondMateSet()) {
      reads.push_back(read_pair.second_mate_ptr_);
    }
  }
}

void ReadPairs::Clear() {
  read_pairs_.clear();
  num_reads_ = 0;
}

}  // namespace reads
