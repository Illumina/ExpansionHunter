//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Peter Krusche <pkrusche@illumina.com>

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

#include "graphalign/GraphAlignment.hh"

#include <memory>

namespace graphtools
{

/**
 * Class to translate between different ways to index nucleotides on the graph
 */
class GraphCoordinates
{
public:
    explicit GraphCoordinates(graphtools::Graph const* g);
    virtual ~GraphCoordinates();
    GraphCoordinates(GraphCoordinates const&) = delete;
    GraphCoordinates& operator=(GraphCoordinates const&) = delete;

    GraphCoordinates(GraphCoordinates&& rhs) noexcept;
    GraphCoordinates& operator=(GraphCoordinates&& rhs) noexcept;

    /**
     * Get a "canonical" / linearized position -- every base on the graph has such a position
     * Positions within the same node are guaranteed to be consecutive
     * @param node node name
     * @param offset offset relative to start of node
     * @return canonical position
     */
    uint64_t canonicalPos(std::string const& node, uint64_t offset = 0) const;

    /**
     * Calculated canonical start and end positions for a graph mapping
     * @param mapping
     * @return start and end
     */
    std::pair<uint64_t, uint64_t> canonicalStartAndEnd(graphtools::Path const& mapping) const;

    /**
     * Reverse lookup : get node and offset from a canonical pos
     * @param canonical_pos canonical position
     * @param node output variable for node name
     * @param offset output variable for offset
     */
    void nodeAndOffset(uint64_t canonical_pos, std::string& node, uint64_t& offset) const;

    /**
     * Calculate the minimum distance in bp between two canonical positions
     * @param pos1 start pos
     * @param pos2 end pos
     * @return basepairs between pos1 and pos2
     */
    uint64_t distance(uint64_t pos1, uint64_t pos2) const;

    /**
     * Return the node id for a node name
     * @param node_name name of node
     * @return node id for the node
     */
    NodeId nodeId(const std::string& node_name) const;

    /**
     * @return the graph for these coordinates
     */
    Graph const& getGraph() const;

private:
    struct GraphCoordinatesImpl;
    std::unique_ptr<GraphCoordinatesImpl> _impl;
};
}
