// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "graphcore/Path.hh"

#include <algorithm>
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

namespace graphtools
{
struct Path::Impl
{
    Impl(
        const Graph* new_graph_raw_ptr, int32_t new_start_position, const vector<NodeId>& new_nodes,
        int32_t new_end_position)
        : graph_raw_ptr(new_graph_raw_ptr)
        , start_position(new_start_position)
        , end_position(new_end_position)
        , nodes(new_nodes)
    {
    }
    bool isValid() const;
    bool isNodePositionValid(NodeId node_id, int32_t position) const;
    bool arePositionsOrdered() const;
    bool isPathEmpty() const;
    bool isFirstNodePosValid() const;
    bool isLastNodePosValid() const;
    bool isPathConnected() const;
    string encode() const;
    void assertThatIndexIsValid(int32_t node_index) const
    {
        if (node_index < 0 || node_index >= (signed)nodes.size())
        {
            const string msg = "Node index " + to_string(node_index) + "is out of bounds for path " + encode();
            throw std::logic_error(msg);
        }
    }

    bool operator==(const Impl& other) const
    {
        return (graph_raw_ptr == other.graph_raw_ptr) && (start_position == other.start_position)
            && (end_position == other.end_position) && (nodes == other.nodes);
    }

    const Graph* const graph_raw_ptr;
    int32_t start_position;
    int32_t end_position;
    vector<NodeId> nodes;
};

bool Path::Impl::isValid() const
{
    return !isPathEmpty() && isFirstNodePosValid() && isLastNodePosValid() && arePositionsOrdered()
        && isPathConnected();
}

bool Path::Impl::arePositionsOrdered() const { return nodes.size() != 1 || start_position <= end_position; }

bool Path::Impl::isNodePositionValid(NodeId node_id, int32_t position) const
{
    if (position < 0)
    {
        return false;
    }
    const string& node_seq = graph_raw_ptr->nodeSeq(node_id);
    return (unsigned)position <= node_seq.length();
}

bool Path::Impl::isPathEmpty() const { return nodes.empty(); }

bool Path::Impl::isFirstNodePosValid() const
{
    const NodeId first_node_id = nodes.front();
    return isNodePositionValid(first_node_id, start_position);
}

bool Path::Impl::isLastNodePosValid() const
{
    const NodeId last_node_id = nodes.back();
    return isNodePositionValid(last_node_id, end_position);
}

bool Path::Impl::isPathConnected() const
{
    vector<NodeId>::const_iterator start_iter;
    vector<NodeId>::const_iterator end_iter;
    for (start_iter = nodes.begin(); start_iter != std::prev(nodes.end()); ++start_iter)
    {
        end_iter = std::next(start_iter);
        if (!graph_raw_ptr->hasEdge(*start_iter, *end_iter))
        {
            return false;
        }
    }
    return true;
}

string Path::Impl::encode() const
{
    string path_encoding;

    size_t node_index = 0;
    const size_t last_index = nodes.size() - 1;
    for (NodeId node_id : nodes)
    {
        const string node_name = to_string(node_id);
        string node_encoding;
        if (node_index == 0) // Encoding first node.
        {
            node_encoding = "(" + node_name + "@" + to_string(start_position) + ")";
        }
        if (node_index == last_index) // Encoding last node.
        {
            node_encoding += "-(" + node_name + "@" + to_string(end_position) + ")";
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

string Path::encode() const { return pimpl_->encode(); }

Path::Path(const Graph* graph_raw_ptr, int32_t start_position, const vector<NodeId>& nodes, int32_t end_position)
    : pimpl_(new Impl(graph_raw_ptr, start_position, nodes, end_position))
{
    if (!pimpl_->isValid())
    {
        throw std::logic_error("Cannot create invalid path");
    }
}

Path::~Path() = default;

Path::Path(const Path& other)
    : pimpl_(new Impl(*other.pimpl_))
{
}

Path::Path(Path&& other) noexcept
    : pimpl_(std::move(other.pimpl_))
{
}

Path& Path::operator=(const Path& other)
{
    if (this != &other)
    {
        pimpl_.reset(new Impl(*other.pimpl_));
    }
    return *this;
}

Path::const_iterator Path::begin() const { return pimpl_->nodes.begin(); }
Path::const_iterator Path::end() const { return pimpl_->nodes.end(); }

Path& Path::operator=(Path&& other) noexcept
{
    pimpl_ = std::move(other.pimpl_);
    return *this;
}

int32_t Path::startPosition() const { return pimpl_->start_position; }
int32_t Path::endPosition() const { return pimpl_->end_position; }
const Graph* Path::graphRawPtr() const { return pimpl_->graph_raw_ptr; }

vector<NodeId> const& Path::nodeIds() const { return pimpl_->nodes; }

size_t Path::numNodes() const { return pimpl_->nodes.size(); }

NodeId Path::getNodeIdByIndex(size_t node_index) const { return pimpl_->nodes[node_index]; }

bool Path::checkOverlapWithNode(NodeId node_id) const
{
    const vector<NodeId>& nodes = pimpl_->nodes;
    return std::find(nodes.begin(), nodes.end(), node_id) != nodes.end();
}

int32_t Path::getStartPositionOnNodeByIndex(size_t node_index) const
{
    pimpl_->assertThatIndexIsValid(static_cast<int32_t>(node_index));

    if (node_index == 0)
    {
        return startPosition();
    }

    return 0;
}

int32_t Path::getEndPositionOnNodeByIndex(size_t node_index) const
{
    pimpl_->assertThatIndexIsValid(static_cast<int32_t>(node_index));

    if (node_index == numNodes() - 1)
    {
        return endPosition();
    }

    const int32_t node_id = pimpl_->nodes[node_index];
    const size_t node_length = graphRawPtr()->nodeSeq(static_cast<NodeId>(node_id)).length();

    return node_length;
}

size_t Path::getNodeOverlapLengthByIndex(size_t node_index) const
{
    pimpl_->assertThatIndexIsValid(static_cast<int32_t>(node_index));
    const int32_t node_id = pimpl_->nodes[node_index];
    const size_t node_length = graphRawPtr()->nodeSeq(static_cast<NodeId>(node_id)).length();
    auto length_on_node = (int32_t)node_length; // This is the length of all intermediate nodes.

    const bool is_first_node = node_index == 0;
    const bool is_last_node = node_index + 1 == numNodes();

    if (is_first_node && is_last_node)
    {
        length_on_node = pimpl_->end_position - pimpl_->start_position;
    }
    else if (is_first_node)
    {
        length_on_node = static_cast<int32_t>(node_length - pimpl_->start_position);
    }
    else if (is_last_node)
    {
        length_on_node = pimpl_->end_position;
    }

    return static_cast<size_t>(length_on_node);
}

int32_t Path::getDistanceFromPathStart(NodeId node, int32_t offset) const
{
    size_t n = 0;
    int32_t distance = 0;
    bool found = false;
    while (n < numNodes())
    {
        const auto node_id = pimpl_->nodes[n];
        const int32_t node_start = n == 0 ? pimpl_->start_position : 0;
        const int32_t node_end
            = n == numNodes() - 1 ? pimpl_->end_position : (int32_t)pimpl_->graph_raw_ptr->nodeSeq(node_id).size() - 1;

        if (node_id == node && offset >= node_start && offset <= node_end)
        {
            distance += offset - node_start;
            found = true;
            break;
        }

        distance += node_end - node_start + 1;
        n++;
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
    for (int32_t node_index = 0; node_index != (signed)pimpl_->nodes.size(); ++node_index)
    {
        path_length += getNodeOverlapLengthByIndex(static_cast<size_t>(node_index));
    }

    return path_length;
}

string Path::getNodeSeq(size_t node_index) const
{
    auto node_id = static_cast<int32_t>(pimpl_->nodes[node_index]);
    const string& sequence = pimpl_->graph_raw_ptr->nodeSeq(static_cast<NodeId>(node_id));

    if (node_index == 0)
    {
        const size_t node_overlap_len = getNodeOverlapLengthByIndex(node_index);
        return sequence.substr(static_cast<unsigned long>(pimpl_->start_position), node_overlap_len);
    }
    else if ((size_t)node_index == pimpl_->nodes.size() - 1)
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
    size_t node_index = 0;
    for (NodeId node_id : pimpl_->nodes)
    {
        string node_seq = pimpl_->graph_raw_ptr->nodeSeq(node_id);
        if (node_index == 0)
        {
            node_seq = node_seq.substr(static_cast<unsigned long>(pimpl_->start_position));
        }

        if (node_index == pimpl_->nodes.size() - 1)
        {
            const int32_t end_node_start = pimpl_->nodes.size() == 1 ? pimpl_->start_position : 0;
            const int32_t segment_len = pimpl_->end_position - end_node_start;
            node_seq = node_seq.substr(0, (unsigned long)segment_len);
        }

        path_seq += node_seq;
        ++node_index;
    }
    return path_seq;
}

bool Path::operator==(const Path& other) const { return *pimpl_ == *other.pimpl_; }

ostream& operator<<(ostream& os, const Path& path) { return os << path.encode(); }

void Path::shiftStartAlongNode(int32_t shift_len)
{
    pimpl_->start_position -= shift_len;
    if (!pimpl_->isValid())
    {
        throw std::logic_error("Cannot move start by " + to_string(shift_len));
    }
}

void Path::shiftEndAlongNode(int32_t shift_len)
{
    pimpl_->end_position += shift_len;

    if (!pimpl_->isValid())
    {
        throw std::logic_error("Cannot move end by " + to_string(shift_len));
    }
}

void Path::extendStartToNode(NodeId node_id)
{
    pimpl_->nodes.insert(pimpl_->nodes.begin(), node_id);
    const auto new_node_seq_len = static_cast<int32_t>(pimpl_->graph_raw_ptr->nodeSeq(node_id).length());
    pimpl_->start_position = new_node_seq_len;

    if (!pimpl_->isValid())
    {
        throw std::logic_error("Cannot extend to node " + to_string(node_id));
    }
}

void Path::extendStartToIncludeNode(NodeId node_id)
{
    pimpl_->nodes.insert(pimpl_->nodes.begin(), node_id);
    pimpl_->start_position = 0;

    if (!pimpl_->isValid())
    {
        throw std::logic_error("Cannot extend to node " + to_string(node_id));
    }
}

void Path::removeStartNode()
{
    pimpl_->nodes.erase(pimpl_->nodes.begin());
    pimpl_->start_position = 0;

    if (!pimpl_->isValid())
    {
        throw std::logic_error("Cannot remove start node of " + encode());
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

void Path::extendEndToNode(NodeId node_id)
{
    pimpl_->nodes.push_back(node_id);
    pimpl_->end_position = 0;

    if (!pimpl_->isValid())
    {
        throw std::logic_error("Cannot extend right to node " + to_string(node_id));
    }
}

void Path::extendEndToIncludeNode(NodeId node_id)
{
    pimpl_->nodes.push_back(node_id);
    const auto new_node_seq_len = static_cast<int32_t>(pimpl_->graph_raw_ptr->nodeSeq(node_id).length());
    pimpl_->end_position = new_node_seq_len;

    if (!pimpl_->isValid())
    {
        throw std::logic_error("Cannot extend right to node " + to_string(node_id));
    }
}

void Path::removeEndNode()
{
    pimpl_->nodes.erase(pimpl_->nodes.end() - 1);
    NodeId new_last_node_id = pimpl_->nodes.back();
    auto new_last_node_len = static_cast<int32_t>(pimpl_->graph_raw_ptr->nodeSeq(new_last_node_id).length());
    pimpl_->end_position = new_last_node_len;
    if (!pimpl_->isValid())
    {
        throw std::logic_error("Cannot remove end node of  " + encode());
    }
}

void Path::shrinkStartBy(int32_t shrink_len)
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

void Path::shrinkEndBy(int32_t shrink_len)
{
    const int32_t node_len_left = pimpl_->end_position;

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

void Path::shrinkBy(int32_t start_shrink_len, int32_t end_shrink_len)
{
    shrinkStartBy(start_shrink_len);
    shrinkEndBy(end_shrink_len);
}

bool Path::operator<(const Path& other) const
{
    if (pimpl_->start_position != other.pimpl_->start_position)
    {
        return pimpl_->start_position < other.pimpl_->start_position;
    }

    if (pimpl_->nodes != other.pimpl_->nodes)
    {
        return pimpl_->nodes < other.pimpl_->nodes;
    }

    return pimpl_->end_position < other.pimpl_->end_position;
}
}
