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

#include "graphs/graph.h"

#include <stdexcept>

using std::set;
using std::string;
using std::to_string;

const string& Graph::NodeSeq(int32_t node_id) const {
  AssertNodeExists(node_id);
  return nodes_[node_id].seq;
}

void Graph::SetNodeSeq(int32_t node_id, const string& seq) {
  AssertNodeExists(node_id);
  nodes_[node_id].seq = seq;
}

void Graph::AddEdge(int32_t source_node_id, int32_t sink_node_id) {
  AssertNodeExists(source_node_id);
  AssertNodeExists(sink_node_id);
  if (source_node_id >= sink_node_id) {
    const string edge =
        "(" + to_string(source_node_id) + "," + to_string(sink_node_id) + ")";
    throw std::logic_error("Edge " + edge + " breaks topological order");
  }

  adjacency_list_[source_node_id].emplace(sink_node_id);
  reverse_adjacency_list_[sink_node_id].emplace(source_node_id);
}

bool Graph::HasEdge(int32_t source_node_id, int32_t sink_node_id) const {
  AssertNodeExists(source_node_id);
  AssertNodeExists(sink_node_id);

  const std::set<int32_t>& successors = adjacency_list_[source_node_id];
  return successors.find(sink_node_id) != successors.end();
}

void Graph::AssertNodeExists(int32_t node_id) const {
  if (node_id < 0 || node_id >= nodes_.size()) {
    throw std::logic_error("Node " + to_string(node_id) + " does not exist");
  }
}

void Graph::Init(int32_t num_nodes) {
  nodes_.resize(num_nodes);
  adjacency_list_.resize(num_nodes);
  reverse_adjacency_list_.resize(num_nodes);
}

void Graph::ClearNodesAndEdges() {
  nodes_.clear();
  adjacency_list_.clear();
  reverse_adjacency_list_.clear();
}