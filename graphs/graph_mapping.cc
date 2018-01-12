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

#include <sstream>
#include <stdexcept>
#include <typeindex>

#include "graphs/graph_mapping.h"

using std::list;
using std::map;
using std::string;
using std::to_string;

void GraphMapping::AssertValidity() const {
  const int32_t mapping_reference_start = mappings_.front().ReferenceStart();

  const Mapping& last_mapping = mappings_.back();
  const int32_t mapping_reference_end =
      last_mapping.ReferenceSpan() + last_mapping.ReferenceStart() - 1;

  if (mapping_reference_start != path_.StartPosition() ||
      mapping_reference_end != path_.EndPosition()) {
    std::ostringstream graph_mapping_encoding;
    graph_mapping_encoding << *this;
    throw std::logic_error("Path " + path_.Encode() +
                           " is not compatible with graph mapping " +
                           graph_mapping_encoding.str());
  }
}

string GraphMapping::Query() const {
  string query;
  for (const auto& mapping : mappings_) {
    query += mapping.Query();
  }
  return query;
}

string GraphMapping::Reference() const {
  string reference;
  for (const auto& mapping : mappings_) {
    reference += mapping.Reference();
  }
  return reference;
}

int32_t GraphMapping::QuerySpan() const {
  int32_t query_span = 0;
  for (const auto& mapping : mappings_) {
    query_span += mapping.QuerySpan();
  }
  return query_span;
}

int32_t GraphMapping::ReferenceSpan() const {
  int32_t reference_span = 0;
  for (const auto& mapping : mappings_) {
    reference_span += mapping.ReferenceSpan();
  }
  return reference_span;
}

int32_t GraphMapping::NumMatches() const {
  int32_t num_matches = 0;
  for (const auto& mapping : mappings_) {
    num_matches += mapping.NumMatched();
  }
  return num_matches;
}

bool GraphMapping::OverlapsNode(int32_t node_id) const {
  return path_.OverlapsNode(node_id);
}

list<int32_t> GraphMapping::GetIndexesOfNode(int32_t node_id) const {
  list<int32_t> indexes;
  const int32_t num_mappings = static_cast<int32_t>(mappings_.size());
  for (int32_t node_index = 0; node_index != num_mappings; ++node_index) {
    const int32_t cur_node_id = path_.GetNodeIdByIndex(node_index);
    if (cur_node_id == node_id) {
      indexes.push_back(node_index);
    }
  }
  return indexes;
}

string GraphMapping::GetCigarString() const {
  string graph_cigar;
  for (int32_t index = 0; index != size(); ++index) {
    const int32_t node_id = path_.GetNodeIdByIndex(index);
    graph_cigar += std::to_string(node_id);
    const Mapping& mapping = mappings_[index];
    graph_cigar += "[" + mapping.GetCigarString() + "]";
  }
  return graph_cigar;
}

std::ostream& operator<<(std::ostream& out, const GraphMapping& graph_mapping) {
  for (int32_t node_index = 0; node_index != graph_mapping.size();
       ++node_index) {
    const int32_t node_id = graph_mapping.GetNodeIdByIndex(node_index);
    out << node_id << "[" << graph_mapping[node_index] << "]";
  }
  return out;
}
