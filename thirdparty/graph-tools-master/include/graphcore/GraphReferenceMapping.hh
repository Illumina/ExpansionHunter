//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>

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

#include <boost/optional.hpp>
#include <memory>
#include <set>

#include "Graph.hh"
#include "Path.hh"

namespace graphtools
{
typedef std::string ContigId; // Identifier of a contig (chromosome) within a reference
typedef int32_t Position; // 0-based position in a reference sequence

/**
 * Defines an interval on a (genomic) reference sequence
 */
class ReferenceInterval
{
public:
    ReferenceInterval(ContigId contig, Position start, Position end);
    /**
     * Create a 0-length interval with start=stop=pos.
     * Represents the position right before base 'pos' (0-based)
     */
    static ReferenceInterval makePosition(ContigId const& contig, Position pos);
    /**
     * Create a region by parsing it from a region string
     * @param region string: <chrName>:<start>-<stop>. 0-based Half-open interval
     * @return ReferenceInterval matching the region
     * @throws If not a valid region string
     */
    static ReferenceInterval parseRegion(std::string const& regionString);
    friend bool operator<(ReferenceInterval const&, ReferenceInterval const&);
    friend bool operator==(ReferenceInterval const&, ReferenceInterval const&);
    friend std::ostream& operator<<(std::ostream&, ReferenceInterval const&);

    /**
     * Length (number of bases covered) of the interval
     */
    int32_t length() const;

    // Reference sequence (chromosome) name
    ContigId const contig;
    // Start 0-based closed (i.e. included)
    Position const start;
    // End 0-based open (i.e. excluded)
    Position const end;
};

/**
 * Map a node to a single piece of reference sequence
 * Very simple 1-1 mapping for now
 */
class NodeReferenceMapping
{
public:
    /**
     * Create a mapping from a node to a reference interval
     * The reference interval must have the same length as the sequence of the node
     */
    NodeReferenceMapping(Graph const&, NodeId, ReferenceInterval const&);

    /**
     * Map a position within a node to a reference position using the NodeReferenceMapping
     * @param offset Position within the node. Must be < node length
     * @return 0-length interval at mapped position
     */
    ReferenceInterval map(int32_t offset) const;

private:
    int32_t const nodeLength_;
    ReferenceInterval const ref_;
};

/**
 * Mapping of (a subset of) nodes in a graph to a reference sequence
 * At most one mapping per node for now
 */
class GraphReferenceMapping
{
public:
    explicit GraphReferenceMapping(Graph const*);

    /**
     * Map a node to a reference interval
     * ReferenceInterval must have the same length as node sequence
     */
    void addMapping(NodeId, ReferenceInterval const&);
    /**
     * Map a position within a node to a reference position
     * @param node in the graph this GraphReferenceMapping is based on
     * @param offset Position within the node. Must be < node length
     * @return 0-length interval at mapped position if the node is mapped, nothing otherwise
     */
    boost::optional<ReferenceInterval> map(NodeId node, int32_t offset) const;
    /**
     * Map the first mappable position in a path to a reference position
     * I.e. the start position of the path in the first node that has a reference mapping
     * @return 0-length interval at mapped position if the node is mapped, nothing otherwise
     */
    boost::optional<ReferenceInterval> map(Path const&) const;

private:
    Graph const* const graph_;
    std::unordered_map<NodeId, NodeReferenceMapping> mappings_;
};
}
