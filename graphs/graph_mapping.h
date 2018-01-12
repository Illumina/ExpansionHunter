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
#include "graphs/path.h"

class GraphMapping {
 public:
  typedef size_t size_type;
  typedef std::vector<Mapping> NodeMappings;
  typedef NodeMappings::const_iterator const_iterator;
  GraphMapping(const GraphPath& path, const std::vector<Mapping>& mappings)
      : path_(path), mappings_(mappings) {
    AssertValidity();
  }
  std::string Query() const;
  std::string Reference() const;
  int32_t QuerySpan() const;
  int32_t ReferenceSpan() const;
  int32_t NumMatches() const;
  bool OverlapsNode(int32_t node_id) const;
  int32_t GetNodeIdByIndex(int32_t node_index) const {
    return path_.GetNodeIdByIndex(node_index);
  }
  std::list<int32_t> GetIndexesOfNode(int32_t node_id) const;
  const_iterator begin() const { return mappings_.begin(); }
  const_iterator end() const { return mappings_.end(); }
  const Mapping& front() const { return mappings_.front(); }
  const Mapping& back() const { return mappings_.back(); }
  size_type size() const { return mappings_.size(); }
  const Mapping& operator[](size_t index) const { return mappings_[index]; }
  bool operator==(const GraphMapping& other) const {
    return path_ == other.path_ && mappings_ == other.mappings_;
  }
  std::string GetCigarString() const;

  friend std::ostream& operator<<(std::ostream& os,
                                  const GraphMapping& graph_mapping);

 private:
  void AssertValidity() const;
  std::vector<Mapping> mappings_;
  GraphPath path_;
};

std::ostream& operator<<(std::ostream& os, const GraphMapping& graph_mapping);

typedef std::unique_ptr<GraphMapping> GraphMappingUPtr;
