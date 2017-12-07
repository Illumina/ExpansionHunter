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
#include <set>
#include <string>
#include <vector>

struct Node {
  std::string name;
  std::string seq;
};

typedef std::vector<std::set<int32_t>> AdjacencyList;

class Graph {
 public:
  Graph(int32_t num_nodes = 0) {
    nodes_.resize(num_nodes);
    adjacency_list_.resize(num_nodes);
    reverse_adjacency_list_.resize(num_nodes);
  }
  int32_t NumNodes() const { return nodes_.size(); }
  void AddEdge(int32_t source_node_id, int32_t sink_node_id);
  bool HasEdge(int32_t source_node_id, int32_t sink_node_id) const;
  const std::set<int32_t>& Successors(int32_t node_id) const {
    AssertNodeExists(node_id);
    return adjacency_list_[node_id];
  }

  const std::set<int32_t>& Predecessors(int32_t node_id) const {
    AssertNodeExists(node_id);
    return reverse_adjacency_list_[node_id];
  }

  void Reset(int32_t num_nodes) {
    nodes_.clear();
    adjacency_list_.clear();
    reverse_adjacency_list_.clear();

    nodes_.resize(num_nodes);
    adjacency_list_.resize(num_nodes);
    reverse_adjacency_list_.resize(num_nodes);
  }

 private:
  void AssertNodeExists(int32_t node_id) const;
  std::vector<Node> nodes_;
  AdjacencyList adjacency_list_;
  AdjacencyList reverse_adjacency_list_;
};