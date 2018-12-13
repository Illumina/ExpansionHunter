//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Peter Krusche <pkrusche@illumina.com>

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
