//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include "graphcore/Path.hh"

#include <algorithm>
#include <iostream>
#include <stdexcept>

using std::logic_error;
using std::ostream;
using std::string;
using std::to_string;
using std::vector;

namespace graphtools
{

void Path::assertValidity() const
{
    assertNonEmpty();
    assertFirstNodePosValid();
    assertLastNodePosValid();
    assertPositionsOrdered();
    assertConnected();
}

void Path::assertPositionsOrdered() const
{
    const bool positionsOrdered = nodes_.size() != 1 || start_position_ <= end_position_;
    if (!positionsOrdered)
    {
        throw logic_error("Positions are not ordered");
    }
}

bool Path::isNodePositionValid(const NodeId node_id, const int32_t position) const
{
    if (position < 0)
    {
        return false;
    }
    const string& node_seq = graph_raw_ptr_->nodeSeq(node_id);
    return (unsigned)position <= node_seq.length();
}

void Path::assertNonEmpty() const
{
    if (nodes_.empty())
    {
        throw logic_error("Path is empty");
    }
}

void Path::assertFirstNodePosValid() const
{
    const NodeId first_node_id = nodes_.front();
    if (!isNodePositionValid(first_node_id, start_position_))
    {
        throw logic_error("Position of first node is invalid");
    }
}

void Path::assertLastNodePosValid() const
{
    const NodeId last_node_id = nodes_.back();
    if (!isNodePositionValid(last_node_id, end_position_))
    {
        throw logic_error("Position of last node is invalid");
    }
}

void Path::assertConnected() const
{
    const unsigned nodeCount(nodes_.size());
    for (unsigned nodeIndex(0); (nodeIndex + 1) < nodeCount; ++nodeIndex)
    {
        if (!graph_raw_ptr_->hasEdge(nodes_[nodeIndex], nodes_[nodeIndex + 1]))
        {
            throw logic_error("Path is not connected");
        }
    }
}

void Path::assertThatIndexIsValid(const int32_t node_index) const
{
    if (node_index < 0 || node_index >= (signed)nodes_.size())
    {
        const string msg = "Node index " + to_string(node_index) + "is out of bounds for path " + encode();
        throw std::logic_error(msg);
    }
}

string Path::encode() const
{
    string path_encoding;

    size_t node_index = 0;
    const size_t last_index = nodes_.size() - 1;
    for (NodeId node_id : nodes_)
    {
        const string node_name = to_string(node_id);
        string node_encoding;
        if (node_index == 0) // Encoding first node.
        {
            node_encoding = "(" + node_name + "@" + to_string(start_position_) + ")";
        }
        if (node_index == last_index) // Encoding last node.
        {
            node_encoding += "-(" + node_name + "@" + to_string(end_position_) + ")";
        }
        if (node_index != 0 && node_index != last_index) // Encoding intermediate node.
        {
            node_encoding = "-(" + node_name + ")";
        }
        path_encoding += node_encoding;
        ++node_index;
    }

    return path_encoding;
}

Path::Path(const Graph* graph_raw_ptr, int32_t start_position, const vector<NodeId>& nodes, int32_t end_position)
    : graph_raw_ptr_(graph_raw_ptr)
    , start_position_(start_position)
    , end_position_(end_position)
    , nodes_(nodes)
{
    try
    {
        assertValidity();
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to create path " + encode() + ": " + e.what());
    }
}

bool Path::checkOverlapWithNode(const NodeId node_id) const
{
    return std::find(nodes_.begin(), nodes_.end(), node_id) != nodes_.end();
}

int32_t Path::getStartPositionOnNodeByIndex(size_t node_index) const
{
    assertThatIndexIsValid(static_cast<int32_t>(node_index));

    if (node_index == 0)
    {
        return startPosition();
    }

    return 0;
}

int32_t Path::getEndPositionOnNodeByIndex(size_t node_index) const
{
    assertThatIndexIsValid(static_cast<int32_t>(node_index));

    if (node_index == (numNodes() - 1))
    {
        return endPosition();
    }

    return graphRawPtr()->nodeSeq(nodes_[node_index]).length();
}

size_t Path::getNodeOverlapLengthByIndex(size_t node_index) const
{
    assertThatIndexIsValid(static_cast<int32_t>(node_index));
    const size_t node_length = graphRawPtr()->nodeSeq(nodes_[node_index]).length();
    auto length_on_node = (int32_t)node_length; // This is the length of all intermediate nodes.

    const bool is_first_node = node_index == 0;
    const bool is_last_node = node_index + 1 == numNodes();

    if (is_first_node && is_last_node)
    {
        length_on_node = end_position_ - start_position_;
    }
    else if (is_first_node)
    {
        length_on_node = static_cast<int32_t>(node_length - start_position_);
    }
    else if (is_last_node)
    {
        length_on_node = end_position_;
    }

    return static_cast<size_t>(length_on_node);
}

int32_t Path::getDistanceFromPathStart(const NodeId node, const int32_t offset) const
{
    int32_t distance = 0;
    bool found = false;
    const size_t nodeCount(numNodes());
    for (size_t nodeIndex(0); nodeIndex < nodeCount; ++nodeIndex)
    {
        const auto node_id = nodes_[nodeIndex];
        const int32_t node_start = nodeIndex == 0 ? start_position_ : 0;
        const int32_t node_end
            = (nodeIndex + 1) == nodeCount ? end_position_ : (int32_t)graph_raw_ptr_->nodeSeq(node_id).size() - 1;

        if (node_id == node && offset >= node_start && offset <= node_end)
        {
            distance += offset - node_start;
            found = true;
            break;
        }

        distance += node_end - node_start + 1;
    }
    if (!found)
    {
        throw std::logic_error(std::to_string(node) + "@" + std::to_string(offset) + " is not on path " + encode());
    }
    return distance;
}

size_t Path::length() const
{
    size_t path_length = 0;
    const size_t nodeCount(numNodes());
    for (size_t node_index = 0; node_index != nodeCount; ++node_index)
    {
        path_length += getNodeOverlapLengthByIndex(node_index);
    }

    return path_length;
}

string Path::getNodeSeq(const size_t node_index) const
{
    const string& sequence = graph_raw_ptr_->nodeSeq(nodes_[node_index]);

    if (node_index == 0)
    {
        const size_t node_overlap_len = getNodeOverlapLengthByIndex(node_index);
        return sequence.substr(static_cast<unsigned long>(start_position_), node_overlap_len);
    }
    else if ((node_index + 1) == numNodes())
    {
        const size_t node_overlap_len = getNodeOverlapLengthByIndex(node_index);
        return sequence.substr(0, node_overlap_len);
    }
    else
    {
        return sequence;
    }
}

string Path::seq() const
{
    string path_seq;
    const size_t nodeCount(numNodes());
    for (size_t node_index = 0; node_index != nodeCount; ++node_index)
    {
        const std::string& node_seq(graph_raw_ptr_->nodeSeq(nodes_[node_index]));
        if ((node_index == 0) or ((node_index + 1) == nodeCount))
        {
            const size_t pos((node_index == 0) ? start_position_ : 0);
            const size_t len(((node_index + 1) == nodeCount) ? (end_position_ - pos) : std::string::npos);
            path_seq += node_seq.substr(pos, len);
        }
        else
        {
            path_seq += node_seq;
        }
    }
    return path_seq;
}

ostream& operator<<(ostream& os, const Path& path) { return os << path.encode(); }

void Path::shiftStartAlongNode(const int32_t shift_len)
{
    start_position_ -= shift_len;

    try
    {
        assertValidity();
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to shift start of " + encode() + " by " + to_string(shift_len) + ": " + e.what());
    }
}

void Path::shiftEndAlongNode(const int32_t shift_len)
{
    end_position_ += shift_len;

    try
    {
        assertValidity();
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to shift end of " + encode() + " by " + to_string(shift_len) + ": " + e.what());
    }
}

void Path::extendStartToNode(const NodeId node_id)
{
    nodes_.insert(nodes_.begin(), node_id);
    const auto new_node_seq_len = static_cast<int32_t>(graph_raw_ptr_->nodeSeq(node_id).length());
    start_position_ = new_node_seq_len;

    try
    {
        assertValidity();
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to extend " + encode() + " to node " + to_string(node_id) + ": " + e.what());
    }
}

void Path::extendStartToIncludeNode(const NodeId node_id)
{
    nodes_.insert(nodes_.begin(), node_id);
    start_position_ = 0;

    try
    {
        assertValidity();
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to extend " + encode() + " to node " + to_string(node_id) + ": " + e.what());
    }
}

void Path::removeStartNode()
{
    nodes_.erase(nodes_.begin());
    start_position_ = 0;

    try
    {
        assertValidity();
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to remove start node of " + encode() + ": " + e.what());
    }
}

void Path::removeZeroLengthStart()
{
    if (numNodes() > 1 && getNodeOverlapLengthByIndex(0) == 0)
    {
        removeStartNode();
    }
}

void Path::removeZeroLengthEnd()
{
    const size_t index_of_last_node = numNodes() - 1;
    if (numNodes() > 1 && getNodeOverlapLengthByIndex(index_of_last_node) == 0)
    {
        removeEndNode();
    }
}

void Path::extendEndToNode(const NodeId node_id)
{
    nodes_.push_back(node_id);
    end_position_ = 0;

    try
    {
        assertValidity();
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to extend " + encode() + " right to node " + to_string(node_id) + ": " + e.what());
    }
}

void Path::extendEndToIncludeNode(const NodeId node_id)
{
    nodes_.push_back(node_id);
    const auto new_node_seq_len = static_cast<int32_t>(graph_raw_ptr_->nodeSeq(node_id).length());
    end_position_ = new_node_seq_len;

    try
    {
        assertValidity();
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to extend " + encode() + " right to node " + to_string(node_id) + ": " + e.what());
    }
}

void Path::removeEndNode()
{
    nodes_.erase(nodes_.end() - 1);
    NodeId new_last_node_id = nodes_.back();
    auto new_last_node_len = static_cast<int32_t>(graph_raw_ptr_->nodeSeq(new_last_node_id).length());
    end_position_ = new_last_node_len;

    try
    {
        assertValidity();
    }
    catch (const std::exception& e)
    {
        throw logic_error("Unable to remove end node of " + encode() + ": " + e.what());
    }
}

void Path::shrinkStartBy(const int32_t shrink_len)
{
    const int32_t node_len_left = getNodeOverlapLengthByIndex(0);

    if (shrink_len <= node_len_left)
    {
        shiftStartAlongNode(-shrink_len);
        removeZeroLengthStart();
    }
    else
    {
        removeStartNode();

        const int32_t leftover_len = shrink_len - node_len_left;
        shrinkStartBy(leftover_len);
    }
}

void Path::shrinkEndBy(const int32_t shrink_len)
{
    const int32_t node_len_left = end_position_;

    if (shrink_len <= node_len_left)
    {
        shiftEndAlongNode(-shrink_len);
        removeZeroLengthEnd();
    }
    else
    {
        removeEndNode();

        const int32_t leftover_len = shrink_len - node_len_left;
        shrinkEndBy(leftover_len);
    }
}

void Path::shrinkBy(const int32_t start_shrink_len, const int32_t end_shrink_len)
{
    shrinkStartBy(start_shrink_len);
    shrinkEndBy(end_shrink_len);
}

bool Path::operator==(const Path& other) const
{
    return (graph_raw_ptr_ == other.graph_raw_ptr_) && (start_position_ == other.start_position_)
        && (end_position_ == other.end_position_) && (nodes_ == other.nodes_);
}

bool Path::operator<(const Path& other) const
{
    if (start_position_ != other.start_position_)
    {
        return start_position_ < other.start_position_;
    }

    if (nodes_ != other.nodes_)
    {
        return nodes_ < other.nodes_;
    }

    return end_position_ < other.end_position_;
}

std::ostream& operator<<(std::ostream& os, const ReversePath& path) { return os << "reverse path of: " << path.path_; }
}
