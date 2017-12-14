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
#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "graphs/graph.h"

// A path in a sequence graph is given by (1) a sequence of nodes and (2)
// start/end position on the first/last node. The start/end positions are
// 0-based and form a closed interval.
class GraphPath {
 public:
  // The constructor does not check if the inputs define a well-formed path;
  // isValid() method can be used to do this.
  GraphPath(std::shared_ptr<Graph> wgraph_ptr, int32_t start_position,
            const std::vector<int32_t>& nodes, int32_t end_position);
  ~GraphPath();
  GraphPath(const GraphPath& other);
  GraphPath(GraphPath&& other) noexcept;
  GraphPath& operator=(const GraphPath& other);
  GraphPath& operator=(GraphPath&& other) noexcept;
  std::vector<int32_t> NodeIds() const;
  size_t NumNodes() const;
  std::string Seq() const;
  std::string SeqOnNodeByIndex(int32_t node_index) const;
  std::shared_ptr<Graph> GraphPtr() const;
  bool IsValid() const;
  bool operator==(const GraphPath& other) const;
  std::string Encode() const;
  int32_t StartPosition() const;
  int32_t EndPosition() const;
  size_t Length() const;
  size_t LengthOnNode(int32_t node_id) const;
  size_t GetOverlapWithNodeByIndex(int32_t node_index) const;
  GraphPath ExtendStartPosition(int32_t extension_len) const;
  GraphPath ExtendEndPosition(int32_t extension_len) const;
  GraphPath ExtendStartNodeTo(int32_t node_id) const;
  GraphPath ExtendEndNodeTo(int32_t node_id) const;
  std::list<GraphPath> ExtendStartBy(int32_t extension_len) const;
  std::list<GraphPath> ExtendEndBy(int32_t extension_len) const;
  // Computes all possible extensions of the path by the specified length in
  // both directions.
  std::list<GraphPath> ExtendBy(int32_t start_extension_len,
                                int32_t end_extension_len) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> pimpl_;
};

std::ostream& operator<<(std::ostream& os, const GraphPath& path);
