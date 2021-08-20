//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Peter Krusche <pkrusche@illumina.com>
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

#include "graphcore/GraphCoordinates.hh"
#include "graphutils/PairHashing.hh"

#include <cassert>

namespace graphtools
{

struct GraphCoordinates::GraphCoordinatesImpl
{
    explicit GraphCoordinatesImpl(Graph const* graph_)
        : graph(graph_)
    {
        uint64_t canonical_offset = 0;
        for (NodeId n_id = 0; n_id < (NodeId)graph_->numNodes(); ++n_id)
        {
            node_name_to_id[graph_->nodeName(n_id)] = n_id;
            auto const& n_name = graph_->nodeName(n_id);
            canonical_offsets[n_name] = canonical_offset;
            node_starts[canonical_offset] = n_name;
            canonical_offset += std::max((size_t)1, graph_->nodeSeq(n_id).size());

            // nodes are sorted in topological order, so we can compute distances as min over all predecessors
            for (NodeId n_source = 0; n_source < (NodeId)graph_->numNodes(); ++n_source)
            {
                // distance = zero in these cases
                if (n_id == n_source || graph_->hasEdge(n_source, n_id))
                {
                    continue;
                }

                size_t min_dist = std::numeric_limits<size_t>::max();
                for (auto pred : graph_->predecessors(n_id))
                {
                    auto pred_distance_it = node_end_to_start_distance.find(std::make_pair(n_source, pred));
                    if (pred_distance_it != node_end_to_start_distance.end())
                    {
                        // minimal distance via that predecessor
                        min_dist = std::min(min_dist, pred_distance_it->second + graph_->nodeSeq(pred).size());
                    }
                    else if (graph_->hasEdge(n_source, pred))
                    {
                        min_dist = std::min(min_dist, graph_->nodeSeq(pred).size());
                    }
                }

                if (min_dist != std::numeric_limits<size_t>::max())
                {
                    node_end_to_start_distance[std::make_pair(n_source, n_id)] = min_dist;
                }
            }
        }
    }

    Graph const* graph;

    std::unordered_map<std::string, uint64_t> canonical_offsets;
    std::map<uint64_t, std::string> node_starts;
    std::map<std::string, NodeId> node_name_to_id;

    std::unordered_map<std::pair<uint64_t, uint64_t>, size_t> node_end_to_start_distance;
};

GraphCoordinates::GraphCoordinates(Graph const* g)
    : _impl(new GraphCoordinatesImpl(g))
{
}
GraphCoordinates::~GraphCoordinates() = default;

GraphCoordinates::GraphCoordinates(GraphCoordinates&& rhs) noexcept
    : _impl(std::move(rhs._impl))
{
}

GraphCoordinates& GraphCoordinates::operator=(GraphCoordinates&& rhs) noexcept
{
    _impl = std::move(rhs._impl);
    return *this;
}

/**
 * Get a "canonical" / linearized position -- every base on the graph has such a position
 * @param node node name
 * @param offset offset relative to start of node
 * @return canonical position
 */
uint64_t GraphCoordinates::canonicalPos(std::string const& node, uint64_t offset) const
{
    auto ioffset = _impl->canonical_offsets.find(node);
    assert(ioffset != _impl->canonical_offsets.end());
    return ioffset->second + offset;
}

/**
 * Calculated canonical start and end positions for a graph mapping
 * @param mapping
 * @return start and end
 */
std::pair<uint64_t, uint64_t> GraphCoordinates::canonicalStartAndEnd(Path const& path) const
{
    std::pair<uint64_t, uint64_t> start_end{ -1, -1 };

    start_end.first
        = canonicalPos(_impl->graph->nodeName(path.nodeIds().front()), static_cast<uint64_t>(path.startPosition()));

    auto end_offset = static_cast<uint64_t>(path.endPosition());
    if (path.numNodes() > 0 && end_offset > 0)
    {
        start_end.second = canonicalPos(_impl->graph->nodeName(path.nodeIds().back()), end_offset);
    }

    if (start_end.first > start_end.second)
    {
        std::swap(start_end.first, start_end.second);
    }

    return start_end;
}

/**
 * Reverse lookup : get node and offset from a canonical pos
 * @param canonical_pos canonical position
 * @param node output variable for node name
 * @param offset output variable for offset
 */
void GraphCoordinates::nodeAndOffset(uint64_t canonical_pos, std::string& node, uint64_t& offset) const
{
    auto lb = _impl->node_starts.lower_bound(canonical_pos);
    if (lb != _impl->node_starts.end())
    {
        if (lb != _impl->node_starts.begin() && canonical_pos < lb->first)
        {
            lb = std::prev(lb);
        }
        node = lb->second;
        offset = canonical_pos - lb->first;
    }
    else
    {
        node = _impl->node_starts.rbegin()->second;
        offset = canonical_pos - _impl->node_starts.rbegin()->first;
    }
}

/**
 * Calculate the minimum distance in bp between two canonical positions
 * @param pos1 start pos
 * @param pos2 end pos
 * @return basepairs between pos1 and pos2
 */
uint64_t GraphCoordinates::distance(uint64_t pos1, uint64_t pos2) const
{
    if (pos1 == pos2)
    {
        return 0;
    }
    if (pos2 < pos1)
    {
        std::swap(pos1, pos2);
    }

    std::string n1, n2;
    uint64_t offset1, offset2;
    nodeAndOffset(pos1, n1, offset1);
    nodeAndOffset(pos2, n2, offset2);

    // on on the same node-> can compute distance directly
    if (n1 == n2)
    {
        return pos2 - pos1;
    }

    const NodeId n1_id = _impl->node_name_to_id[n1];
    const NodeId n2_id = _impl->node_name_to_id[n2];
    const size_t n1_length = _impl->graph->nodeSeq(n1_id).size();

    uint64_t result = std::numeric_limits<uint64_t>::max();
    if (_impl->graph->hasEdge(n1_id, n2_id))
    {
        result = n1_length - offset1 + offset2;
    }
    else
    {
        auto dist_it = _impl->node_end_to_start_distance.find(std::make_pair(n1_id, n2_id));

        if (dist_it != _impl->node_end_to_start_distance.end())
        {
            result = n1_length - offset1 + offset2 + dist_it->second;
        }
    }
    return result;
}

/**
 * Return the node id for a node name
 * @param node_name name of node
 * @return node id for the node
 */
NodeId GraphCoordinates::nodeId(const std::string& node_name) const
{
    assert(_impl->node_name_to_id.count(node_name) != 0);
    return _impl->node_name_to_id[node_name];
}

/**
 * @return the graph for these coordinates
 */
Graph const& GraphCoordinates::getGraph() const { return *(_impl->graph); }
}
