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
                           "  and " + sequence + " have same length");
  }

  vector<string> split_seq;
  const vector<int32_t> path_node_ids = path.NodeIds();

  size_t cur_position = 0;
  for (int32_t node_id : path_node_ids) {
    const size_t length_on_node = path.LengthOnNode(node_id);
    split_seq.push_back(sequence.substr(cur_position, length_on_node));
    cur_position += length_on_node;
  }
  return split_seq;
}
