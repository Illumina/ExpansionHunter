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

#include "graphcore/Graph.hh"

#include <stdexcept>

#include "graphutils/SequenceOperations.hh"

using std::set;
using std::string;
using std::to_string;
using std::vector;

namespace graphtools
{

void Graph::init(size_t numnodes_)
{
    nodes_.resize(numnodes_);
    adjacency_list_.resize(numnodes_);
    reverse_adjacency_list_.resize(numnodes_);
}

void Graph::assertNodeExists(NodeId node_id) const
{
    if (node_id >= nodes_.size())
    {
        throw std::logic_error("Node with id " + to_string(node_id) + " does not exist");
    }
}

void Graph::assertEdgeExists(NodeIdPair node_id_pair) const
{
    if (!hasEdge(node_id_pair.first, node_id_pair.second))
    {
        throw std::logic_error(
            "There is no edge between " + to_string(node_id_pair.first) + " and " + to_string(node_id_pair.second));
    }
}

void assertValidSequence(string const& seq)
{
    if (seq.empty())
    {
        throw std::logic_error("Invalid node sequence " + seq);
    }
}

const string& Graph::nodeName(NodeId node_id) const
{
    assertNodeExists(node_id);
    return nodes_[node_id].name;
}
void Graph::setNodeName(NodeId node_id, const std::string& node_name)
{
    assertNodeExists(node_id);
    nodes_[node_id].name = node_name;
}

const string& Graph::nodeSeq(NodeId node_id) const
{
    assertNodeExists(node_id);
    return nodes_[node_id].sequence;
}

const vector<string>& Graph::nodeSeqExpansion(NodeId node_id) const
{
    assertNodeExists(node_id);
    return nodes_[node_id].sequence_expansion;
}

void Graph::setNodeSeq(NodeId node_id, const string& sequence)
{
    assertNodeExists(node_id);
    assertValidSequence(sequence);
    nodes_[node_id].sequence = sequence;
    if (is_sequence_expansion_required_)
    {
        expandReferenceSequence(sequence, nodes_[node_id].sequence_expansion);
    }
    else
    {
        nodes_[node_id].sequence_expansion = { sequence };
    }
}

void Graph::addEdge(NodeId source_id, NodeId sink_id)
{
    assertNodeExists(source_id);
    assertNodeExists(sink_id);

    const string edge_encoding = "(" + to_string(source_id) + " ," + to_string(sink_id) + ")";
    if (hasEdge(source_id, sink_id))
    {
        throw std::logic_error("Graph already contains edge " + edge_encoding);
    }

    if (source_id > sink_id)
    {
        throw std::logic_error("Edge " + edge_encoding + " breaks topological order");
    }

    NodeIdPair node_id_pair = std::make_pair(source_id, sink_id);
    edge_labels_[node_id_pair];

    adjacency_list_[source_id].emplace(sink_id);
    reverse_adjacency_list_[sink_id].emplace(source_id);
}

bool Graph::hasEdge(NodeId source_id, NodeId sink_id) const
{
    assertNodeExists(source_id);
    assertNodeExists(sink_id);
    NodeIdPair node_id_pair(source_id, sink_id);
    return edge_labels_.find(node_id_pair) != edge_labels_.end();
}

void Graph::addLabelToEdge(NodeId source_id, NodeId sink_id, const std::string& label)
{
    NodeIdPair node_pair = std::make_pair(source_id, sink_id);
    assertEdgeExists(node_pair);
    edge_labels_[node_pair].insert(label);
}

SortedLabels Graph::allLabels() const
{
    SortedLabels labels;
    for (const auto& single_edge_labels_ : edge_labels_)
    {
        labels.insert(single_edge_labels_.second.begin(), single_edge_labels_.second.end());
    }
    return labels;
}

const Labels& Graph::edgeLabels(NodeId source_id, NodeId sink_id) const
{
    NodeIdPair node_id_pair = std::make_pair(source_id, sink_id);
    return edge_labels_.at(node_id_pair);
}

std::set<NodeIdPair> Graph::edgesWithLabel(const std::string& label) const
{
    std::set<NodeIdPair> result;
    for (const auto& pair : edge_labels_)
    {
        if (pair.second.count(label))
        {
            result.insert(pair.first);
        }
    }
    return result;
}

void Graph::eraseLabel(const std::string& label)
{
    for (auto& pair : edge_labels_)
    {
        pair.second.erase(label);
    }
}

const set<NodeId>& Graph::successors(NodeId node_id) const
{
    assertNodeExists(node_id);
    return adjacency_list_[node_id];
}

const std::set<NodeId>& Graph::predecessors(NodeId node_id) const
{
    assertNodeExists(node_id);
    return reverse_adjacency_list_[node_id];
}
}
