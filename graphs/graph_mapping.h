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

#include <initializer_list>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include "graphs/graph.h"
#include "graphs/linear_mapping.h"

struct NodeMapping {
  int32_t node_id;
  Mapping mapping;
  bool operator==(const NodeMapping& other) const {
    return node_id == other.node_id && mapping == other.mapping;
  }
  std::string GetCigarString() const;
};

class GraphMapping {
 public:
  typedef size_t size_type;
  typedef std::vector<NodeMapping> NodeMappings;
  typedef NodeMappings::const_iterator const_iterator;
  GraphMapping() = default;
  GraphMapping(const std::vector<int32_t>& node_ids,
               const std::vector<Mapping>& mappings) {
    for (size_t index = 0; index != node_ids.size(); ++index) {
      NodeMapping node_mapping;
      node_mapping.node_id = node_ids[index];
      node_mapping.mapping = mappings[index];
      node_mappings_.push_back(node_mapping);
    }
  }
  std::string Query() const;
  std::string Reference() const;
  int32_t QuerySpan() const;
  int32_t ReferenceSpan() const;
  int32_t NumMatches() const;
  bool OverlapsNode(int32_t node_id) const;
  std::list<int32_t> GetIndexesOfNode(int32_t node_id) const;
  const_iterator begin() const { return node_mappings_.begin(); }
  const_iterator end() const { return node_mappings_.end(); }
  const NodeMapping& front() const { return node_mappings_.front(); }
  const NodeMapping& back() const { return node_mappings_.back(); }
  size_type size() const { return node_mappings_.size(); }
  const NodeMapping& operator[](size_t index) const {
    return node_mappings_[index];
  }
  bool operator==(const GraphMapping& other) const {
    return node_mappings_ == other.node_mappings_;
  }
  std::string GetCigarString() const;

  friend std::ostream& operator<<(std::ostream& os,
                                  const GraphMapping& graph_mapping);

 private:
  NodeMappings node_mappings_;
};

std::ostream& operator<<(std::ostream& os, const GraphMapping& graph_mapping);

typedef std::unique_ptr<GraphMapping> GraphMappingPtr;
