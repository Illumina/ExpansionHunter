//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>

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

#pragma once

#include <algorithm>
#include <cstdint>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "graphutils/PairHashing.hh"

namespace graphtools
{

using NodeId = uint32_t;
using NodeIdPair = std::pair<NodeId, NodeId>;
using SortedLabels = std::set<std::string>;
using Labels = std::unordered_set<std::string>;
using AdjacencyList = std::vector<std::set<NodeId>>;

struct Node
{
    std::string name;
    std::string sequence;
    std::vector<std::string> sequence_expansion;
};

/**
 * Sequence graph that can hold degenerate nucleotide sequences
 */
class Graph
{
public:
    explicit Graph(size_t num_nodes = 0, std::string const& id = "")
        : graphId(id)
    {
        init(num_nodes);
    }

    virtual ~Graph() = default;

    size_t numNodes() const { return nodes_.size(); }
    size_t numEdges() const { return edge_labels_.size(); }
    const std::string& nodeName(NodeId node_id) const;
    void setNodeName(NodeId node_id, const std::string& node_name);
    const std::string& nodeSeq(NodeId node_id) const;
    void setNodeSeq(NodeId node_id, const std::string& sequence);
    const std::vector<std::string>& nodeSeqExpansion(NodeId node_id) const;
    void addEdge(NodeId source_node_id, NodeId sink_node_id);
    bool hasEdge(NodeId source_node_id, NodeId sink_node_id) const;

    SortedLabels allLabels() const;
    const Labels& edgeLabels(NodeId source_node_id, NodeId sink_node_id) const;
    // All edges in the graph with this label
    std::set<NodeIdPair> edgesWithLabel(const std::string& label) const;
    /**
     * All label to an existing edge
     * @throws if source -> sink is not an edge in the graph
     */
    void addLabelToEdge(NodeId source_node_id, NodeId sink_node_id, const std::string& label);
    // Remove this label from all edges
    void eraseLabel(const std::string& label);

    // this cannot be const if graphs are to be assigned. Currently this happens in unit tests.
    std::string graphId;
    const std::set<NodeId>& successors(NodeId node_id) const;
    const std::set<NodeId>& predecessors(NodeId node_id) const;

private:
    void init(size_t num_nodes);
    void assertNodeExists(NodeId node_id) const;
    void assertEdgeExists(NodeIdPair edge) const;
    std::vector<Node> nodes_;
    std::unordered_map<NodeIdPair, Labels> edge_labels_;
    AdjacencyList adjacency_list_;
    AdjacencyList reverse_adjacency_list_;
};

class ReverseGraph
{
    const Graph& graph_;

public:
    explicit ReverseGraph(const Graph& graph)
        : graph_(graph)
    {
    }

    size_t numNodes() const { return graph_.numNodes(); }
    size_t numEdges() const { return graph_.numEdges(); }
    const std::string& nodeName(NodeId nodeId) const { return graph_.nodeName(nodeId); }
    const std::string nodeSeq(NodeId nodeId) const
    {
        std::string ret = graph_.nodeSeq(nodeId);
        std::reverse(ret.begin(), ret.end());
        return ret;
    }
    bool hasEdge(NodeId sourceNodeId, NodeId sinkNodeId) const { return graph_.hasEdge(sourceNodeId, sinkNodeId); }

    SortedLabels allLabels() const { return graph_.allLabels(); }
    const Labels& edgeLabels(NodeId sourceNodeId, NodeId sinkNodeId) const
    {
        return graph_.edgeLabels(sourceNodeId, sinkNodeId);
    }

    // All edges in the graph with this label
    std::set<NodeIdPair> edgesWithLabel(const std::string& label) const { return graph_.edgesWithLabel(label); }

    const std::set<NodeId>& successors(NodeId nodeId) const { return graph_.predecessors(nodeId); }
    const std::set<NodeId>& predecessors(NodeId nodeId) const { return graph_.successors(nodeId); }
};
}
