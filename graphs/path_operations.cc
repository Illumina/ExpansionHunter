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

#include "graphs/path_operations.h"

using std::list;
using std::string;
using std::vector;

vector<string> SplitByPath(const GraphPath& path, const std::string& sequence) {
  if (path.Length() != sequence.length()) {
    throw std::logic_error("Split operation requires that " + path.Encode() +
                           " and " + sequence + " have same length");
  }

  vector<string> split_seq;
  size_t cur_position = 0;
  for (int32_t node_index = 0; node_index != path.NumNodes(); ++node_index) {
    const size_t length_on_node = path.GetOverlapWithNodeByIndex(node_index);
    split_seq.push_back(sequence.substr(cur_position, length_on_node));
    cur_position += length_on_node;
  }
  return split_seq;
}

list<GraphPath> ComputeRightEndings(const GraphPath& path,
                                    int32_t dist_from_right_end) {
  GraphPath shortened_path = path.ShrinkEndBy(dist_from_right_end);
  const int32_t last_node_index = shortened_path.NumNodes() - 1;
  const int32_t last_node_id = shortened_path.GetNodeIdByIndex(last_node_index);
  const int32_t end_position = shortened_path.EndPosition();
  GraphPath seed_path(path.GraphPtr(), end_position, {last_node_id},
                      end_position);
  return seed_path.ExtendEndBy(dist_from_right_end);
}

list<GraphPath> ComputeLeftEndings(const GraphPath& path,
                                   int32_t dist_from_left_end) {
  GraphPath shortened_path = path.ShrinkStartBy(dist_from_left_end);
  const int32_t first_node_id = shortened_path.GetNodeIdByIndex(0);
  const int32_t start_position = shortened_path.StartPosition();
  GraphPath seed_path(path.GraphPtr(), start_position, {first_node_id},
                      start_position);
  return seed_path.ExtendStartBy(dist_from_left_end);
}
