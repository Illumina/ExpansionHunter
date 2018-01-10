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

#include <stdexcept>
#include <typeindex>

#include "graphs/graph_mapping.h"

using std::list;
using std::map;
using std::string;
using std::to_string;

string NodeMapping::GetCigarString() const {
  string cigar_string = std::to_string(node_id);
  cigar_string += "[" + mapping.GetCigarString() + "]";
  return cigar_string;
}

string GraphMapping::Query() const {
  string query;
  for (const auto& node_mapping : node_mappings_) {
    query += node_mapping.mapping.Query();
  }
  return query;
}

string GraphMapping::Reference() const {
  string reference;
  for (const auto& node_mapping : node_mappings_) {
    reference += node_mapping.mapping.Reference();
  }
  return reference;
}

int32_t GraphMapping::QuerySpan() const {
  int32_t query_span = 0;
  for (const auto& node_mapping : node_mappings_) {
    query_span += node_mapping.mapping.QuerySpan();
  }
  return query_span;
}

int32_t GraphMapping::ReferenceSpan() const {
  int32_t reference_span = 0;
  for (const auto& node_mapping : node_mappings_) {
    reference_span += node_mapping.mapping.ReferenceSpan();
  }
  return reference_span;
}

int32_t GraphMapping::NumMatches() const {
  int32_t num_matches = 0;
  for (const auto& node_mapping : node_mappings_) {
    num_matches += node_mapping.mapping.NumMatched();
  }
  return num_matches;
}

bool GraphMapping::OverlapsNode(int32_t node_id) const {
  for (const auto& node_mapping : node_mappings_) {
    if (node_mapping.node_id == node_id) {
      return true;
    }
  }
  return false;
}

list<int32_t> GraphMapping::GetIndexesOfNode(int32_t node_id) const {
  list<int32_t> indexes;
  for (int32_t node_index = 0; node_index != (int32_t)node_mappings_.size();
       ++node_index) {
    if (node_mappings_[node_index].node_id == node_id) {
      indexes.push_back(node_index);
    }
  }
  return indexes;
}

string GraphMapping::GetCigarString() const {
  string cigar_string;
  for (const auto& node_mapping : node_mappings_) {
    cigar_string += node_mapping.GetCigarString();
  }
  return cigar_string;
}

std::ostream& operator<<(std::ostream& os, const GraphMapping& graph_mapping) {
  for (const NodeMapping& node_mapping : graph_mapping) {
    os << node_mapping.node_id << "[" << node_mapping.mapping << "]";
  }
  return os;
}
