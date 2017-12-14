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

#include "graphs/gapless_aligner.h"

#include <list>
#include <stdexcept>
#include <vector>

#include "graphs/path_operations.h"

using std::list;
using std::string;
using std::vector;

Mapping AlignWithoutGaps(const std::string& query, int32_t ref_start,
                         const std::string& reference) {
  if (reference.length() < ref_start + query.length()) {
    const string msg = "Gapless alignment requires that read " + query +
                       " is shorter than reference " + reference;
    throw std::logic_error(msg);
  }

  if (query.empty() || reference.empty()) {
    throw std::logic_error("Cannot align empty sequences");
  }

  vector<Operation> operations;
  int32_t previous_run_end = 0;
  int32_t run_length = 0;
  char run_operation = '\0';
  for (size_t index = 0; index != query.length(); ++index) {
    char cur_operation = 'X';
    if (query[index] == reference[ref_start + index]) {
      cur_operation = 'M';
    }

    if (cur_operation == run_operation) {
      ++run_length;
    } else {
      if (run_operation != '\0') {
        const string query_piece = query.substr(previous_run_end, run_length);
        const string reference_piece =
            reference.substr(ref_start + previous_run_end, run_length);
        operations.push_back(
            Operation(run_operation, run_length, query_piece, reference_piece));
      }
      previous_run_end += run_length;
      run_length = 1;
      run_operation = cur_operation;
    }
  }

  const string query_piece = query.substr(previous_run_end, run_length);
  const string reference_piece =
      reference.substr(ref_start + previous_run_end, run_length);
  operations.push_back(
      Operation(run_operation, run_length, query_piece, reference_piece));

  return Mapping(ref_start, operations);
}

GraphMapping AlignWithoutGaps(const GraphPath& path, const string& read) {
  vector<string> read_pieces = SplitByPath(path, read);
  vector<Mapping> node_mappings;

  size_t index = 0;
  for (int32_t node_id : path.NodeIds()) {
    std::shared_ptr<Graph> graph_ptr = path.GraphPtr();
    const string node_seq = graph_ptr->NodeSeq(node_id);
    const string read_piece = read_pieces[index];
    const int32_t ref_start = index == 0 ? path.StartPosition() : 0;
    node_mappings.push_back(AlignWithoutGaps(read_piece, ref_start, node_seq));
    ++index;
  }

  return GraphMapping(path.NodeIds(), node_mappings);
}
