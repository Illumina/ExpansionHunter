//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

const string& Graph::nodeSeq(const NodeId node_id) const
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
    expandReferenceSequence(sequence, nodes_[node_id].sequence_expansion);
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
