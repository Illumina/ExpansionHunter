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

#include "graphs/path.h"

#include <cassert>
#include <set>
#include <sstream>
#include <stdexcept>

using std::list;
using std::ostream;
using std::set;
using std::shared_ptr;
using std::string;
using std::to_string;
using std::vector;

struct GraphPath::Impl {
  Impl(shared_ptr<Graph> graph_ptr, int32_t start_position,
       const vector<int32_t>& nodes, int32_t end_position)
      : graph_ptr_(graph_ptr),
        start_position_(start_position),
        end_position_(end_position),
        nodes_(nodes) {}
  bool isNodePositionValid(int32_t node_id, int32_t position) const;
  bool areNodesOrdered() const;
  bool arePositionsOrdered() const;
  bool isPathEmpty() const;
  bool isFirstNodePosValid() const;
  bool isLastNodePosValid() const;
  bool isPathConected() const;

  bool operator==(const Impl& other) const {
    return (graph_ptr_ == other.graph_ptr_) &&
           (start_position_ == other.start_position_) &&
           (end_position_ == other.end_position_) && (nodes_ == other.nodes_);
  }

  shared_ptr<Graph> graph_ptr_;
  int32_t start_position_;
  int32_t end_position_;
  vector<int32_t> nodes_;
};

bool GraphPath::Impl::areNodesOrdered() const {
  int32_t cur_node_id = nodes_.front();
  vector<int32_t>::const_iterator node_id_iter = nodes_.begin();
  ++node_id_iter;  // Assuming the path contains at least one node.
  while (node_id_iter != nodes_.end()) {
    const int32_t next_node_id = *node_id_iter;
    if (cur_node_id >= next_node_id) {
      return false;
    }
    cur_node_id = next_node_id;
    ++node_id_iter;
  }
  return true;
}

bool GraphPath::Impl::arePositionsOrdered() const {
  if (nodes_.size() == 1 && start_position_ > end_position_) {
    return false;
  }
  return true;
}

bool GraphPath::Impl::isNodePositionValid(int32_t node_id,
                                          int32_t position) const {
  if (position < 0) {
    return false;
  }
  const string& node_seq = graph_ptr_->NodeSeq(node_id);
  if ((int32_t)node_seq.length() <= position) {
    return false;
  }
  return true;
}

bool GraphPath::Impl::isPathEmpty() const { return nodes_.empty(); }

bool GraphPath::Impl::isFirstNodePosValid() const {
  const int32_t first_node_id = nodes_.front();
  return isNodePositionValid(first_node_id, start_position_);
}

bool GraphPath::Impl::isLastNodePosValid() const {
  const int32_t last_node_id = nodes_.back();
  return isNodePositionValid(last_node_id, end_position_);
}

bool GraphPath::Impl::isPathConected() const {
  vector<int32_t>::const_iterator start_iter;
  vector<int32_t>::const_iterator end_iter;
  for (start_iter = nodes_.begin(); start_iter != std::prev(nodes_.end());
       ++start_iter) {
    end_iter = std::next(start_iter);
    if (!graph_ptr_->HasEdge(*start_iter, *end_iter)) {
      return false;
    }
  }
  return true;
}

GraphPath::GraphPath(shared_ptr<Graph> graph_ptr, int32_t start_position,
                     const vector<int32_t>& nodes, int32_t end_position)
    : pimpl_(new Impl(graph_ptr, start_position, nodes, end_position)) {}

GraphPath::~GraphPath() = default;

GraphPath::GraphPath(const GraphPath& other)
    : pimpl_(new Impl(*other.pimpl_)) {}

GraphPath::GraphPath(GraphPath&& other) noexcept
    : pimpl_(std::move(other.pimpl_)) {}

GraphPath& GraphPath::operator=(const GraphPath& other) {
  if (this != &other) {
    pimpl_.reset(new Impl(*other.pimpl_));
  }
  return *this;
}

GraphPath& GraphPath::operator=(GraphPath&& other) noexcept {
  pimpl_ = std::move(other.pimpl_);
  return *this;
}

int32_t GraphPath::start_position() const { return pimpl_->start_position_; }
int32_t GraphPath::end_position() const { return pimpl_->end_position_; }
std::shared_ptr<Graph> GraphPath::graph_ptr() const {
  return pimpl_->graph_ptr_;
}

vector<int32_t> GraphPath::node_ids() const { return pimpl_->nodes_; }

size_t GraphPath::num_nodes() const { return pimpl_->nodes_.size(); }

size_t GraphPath::lengthOnNode(int32_t node_id) const {
  const size_t node_length = pimpl_->graph_ptr_->NodeSeq(node_id).length();
  size_t length_on_node =
      node_length;  // This is the length of all intermediate nodes.

  const bool is_first_node = node_id == pimpl_->nodes_.front();
  const bool is_last_node = node_id == pimpl_->nodes_.back();

  if (is_first_node && is_last_node) {
    length_on_node = pimpl_->end_position_ - pimpl_->start_position_ + 1;
  } else if (is_first_node) {
    length_on_node = node_length - pimpl_->start_position_;
  } else if (is_last_node) {
    length_on_node = pimpl_->end_position_ + 1;
  }

  return length_on_node;
}

size_t GraphPath::length() const {
  size_t path_length = 0;
  for (int32_t node_id : pimpl_->nodes_) {
    path_length += lengthOnNode(node_id);
  }

  return path_length;
}

string GraphPath::seqOnNode(int32_t node_index) const {
  uint64_t node_id = (int32_t)pimpl_->nodes_[node_index];
  const string& sequence = pimpl_->graph_ptr_->NodeSeq(node_id);

  if (node_index == 0) {
    return sequence.substr(pimpl_->start_position_, lengthOnNode(node_id));
  } else if ((size_t)node_index == pimpl_->nodes_.size() - 1) {
    return sequence.substr(0, lengthOnNode(node_id));
  } else {
    return sequence;
  }
}

string GraphPath::seq() const {
  string path_seq;
  size_t node_index = 0;
  for (int32_t node_id : pimpl_->nodes_) {
    string node_seq = pimpl_->graph_ptr_->NodeSeq(node_id);
    if (node_index == 0) {
      node_seq = node_seq.substr(pimpl_->start_position_);
    }

    if (node_index == pimpl_->nodes_.size() - 1) {
      const int32_t end_node_start =
          pimpl_->nodes_.size() == 1 ? pimpl_->start_position_ : 0;
      const int32_t segment_len = pimpl_->end_position_ - end_node_start + 1;
      node_seq = node_seq.substr(0, segment_len);
    }

    path_seq += node_seq;
    ++node_index;
  }
  return path_seq;
}

bool GraphPath::isValid() const {
  return (!pimpl_->isPathEmpty() && pimpl_->isFirstNodePosValid() &&
          pimpl_->isLastNodePosValid() && pimpl_->areNodesOrdered() &&
          pimpl_->arePositionsOrdered() && pimpl_->isPathConected());
}

bool GraphPath::operator==(const GraphPath& other) const {
  return *pimpl_ == *other.pimpl_;
}

string GraphPath::encode() const {
  string path_encoding;

  size_t node_index = 0;
  const size_t last_index = pimpl_->nodes_.size() - 1;
  for (int32_t node_id : pimpl_->nodes_) {
    const string node_name = to_string(node_id);
    string node_encoding;
    if (node_index == 0)  // Encoding first node.
    {
      node_encoding =
          "(" + node_name + "@" + to_string(pimpl_->start_position_) + ")";
    }

    if (node_index == last_index)  // Encoding last node.
    {
      node_encoding +=
          "-(" + node_name + "@" + to_string(pimpl_->end_position_) + ")";
    }

    if (node_index != 0 &&
        node_index != last_index)  // Encoding intermediate node.
    {
      node_encoding = "-(" + node_name + ")";
    }
    path_encoding += node_encoding;
    ++node_index;
  }

  return path_encoding;
}

ostream& operator<<(ostream& os, const GraphPath& path) {
  return os << path.encode();
}

GraphPath GraphPath::extendStartPosition(int32_t extension_len) const {
  const int32_t extended_first_node_pos =
      pimpl_->start_position_ - extension_len;
  GraphPath extended_path(pimpl_->graph_ptr_, extended_first_node_pos,
                          pimpl_->nodes_, pimpl_->end_position_);
  if (!extended_path.isValid()) {
    throw std::logic_error("Cannot extend " + encode() + " left by " +
                           to_string(extension_len));
  }
  return extended_path;
}

GraphPath GraphPath::extendEndPosition(int32_t extension_len) const {
  int32_t extended_last_node_pos = pimpl_->end_position_ + extension_len;
  GraphPath extended_path(pimpl_->graph_ptr_, pimpl_->start_position_,
                          pimpl_->nodes_, extended_last_node_pos);
  if (!extended_path.isValid()) {
    throw std::logic_error("Cannot extend " + encode() + " right by " +
                           to_string(extension_len));
  }
  return extended_path;
}

GraphPath GraphPath::extendStartNodeTo(int32_t node_id) const {
  vector<int32_t> extended_nodes = pimpl_->nodes_;
  extended_nodes.insert(extended_nodes.begin(), node_id);
  int32_t new_node_seq_len = pimpl_->graph_ptr_->NodeSeq(node_id).length();
  GraphPath extended_path(pimpl_->graph_ptr_, new_node_seq_len - 1,
                          extended_nodes, pimpl_->end_position_);
  if (!extended_path.isValid()) {
    throw std::logic_error("Cannot extend " + encode() + " left to node %i" +
                           to_string(node_id));
  }

  return extended_path;
}

GraphPath GraphPath::extendEndNodeTo(int32_t node_id) const {
  vector<int32_t> extended_nodes = pimpl_->nodes_;
  extended_nodes.push_back(node_id);
  GraphPath extended_path(pimpl_->graph_ptr_, pimpl_->start_position_,
                          extended_nodes, 0);
  if (!extended_path.isValid()) {
    throw std::logic_error("Cannot extend " + encode() + " right to node " +
                           to_string(node_id));
  }

  return extended_path;
}

list<GraphPath> GraphPath::extendStartBy(int32_t extension_len) const {
  list<GraphPath> extended_paths;

  const int32_t start_node_id = pimpl_->nodes_.front();
  // const int32_t start_node_length =
  // graph_ptr_->node(start_node_id)->sequence().length();

  // Start position gives the maximum extension.
  if (extension_len <= pimpl_->start_position_) {
    extended_paths.push_back(extendStartPosition(extension_len));
  } else {
    const set<int32_t> pred_node_ids =
        pimpl_->graph_ptr_->Predecessors(start_node_id);
    const int32_t leftover_length = extension_len - pimpl_->start_position_ - 1;
    for (int32_t pred_node_id : pred_node_ids) {
      const GraphPath path_with_this_node = extendStartNodeTo(pred_node_id);
      list<GraphPath> extensions_of_path_with_this_node =
          path_with_this_node.extendStartBy(leftover_length);
      extended_paths.splice(extended_paths.end(),
                            extensions_of_path_with_this_node);
    }
  }

  return extended_paths;
}

list<GraphPath> GraphPath::extendEndBy(int32_t extension_len) const {
  list<GraphPath> extended_paths;

  const int32_t end_node_id = pimpl_->nodes_.back();
  const int32_t end_node_length =
      pimpl_->graph_ptr_->NodeSeq(end_node_id).length();
  const int32_t max_extension_at_end_node =
      end_node_length - pimpl_->end_position_ - 1;

  if (extension_len <= max_extension_at_end_node) {
    extended_paths.push_back(extendEndPosition(extension_len));
  } else {
    const set<int32_t> succ_node_ids =
        pimpl_->graph_ptr_->Successors(end_node_id);
    const int32_t leftover_length =
        extension_len - max_extension_at_end_node - 1;
    for (int32_t succ_node_id : succ_node_ids) {
      const GraphPath path_with_this_node = extendEndNodeTo(succ_node_id);
      list<GraphPath> extensions_of_path_with_this_node =
          path_with_this_node.extendEndBy(leftover_length);
      extended_paths.splice(extended_paths.end(),
                            extensions_of_path_with_this_node);
    }
  }

  return extended_paths;
}

list<GraphPath> GraphPath::extendBy(int32_t start_extension_len,
                                    int32_t end_extension_len) const {
  list<GraphPath> extended_paths;
  list<GraphPath> start_extended_paths = extendStartBy(start_extension_len);
  for (const GraphPath& path : start_extended_paths) {
    list<GraphPath> end_extended_paths = path.extendEndBy(end_extension_len);
    extended_paths.splice(extended_paths.end(), end_extended_paths);
  }
  return extended_paths;
}