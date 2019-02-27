//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>
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

#pragma once

#include <cstdint>
#include <iostream>
#include <list>
#include <memory>
#include <string>
#include <vector>

#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

namespace graphtools
{

/**
 * Defines a Path Family (a set of paths) from a set of edges.
 * Briefly, a path is part of a path family F if
 *  - it uses at least one edge in F
 *  - it uses an edge in F into or out-of any node where such an edge is present
 */
class PathFamily
{
public:
    explicit PathFamily(Graph*);
    // Path family from all edges with the given label
    PathFamily(Graph* graph, const std::string& label);
    PathFamily(const PathFamily&);
    PathFamily(PathFamily&&) noexcept;
    ~PathFamily() noexcept;

    bool operator==(const PathFamily&) const;
    bool operator!=(const PathFamily&) const;

    std::unordered_set<NodeIdPair> const& edges() const;
    const Graph& graph() const;

    // Check if a path in contained in (consistent with) F
    bool containsPath(const Path&) const;
    // Check if another path family contains a subset of edges in F
    bool includes(const PathFamily&) const;

    void addEdge(NodeId first, NodeId second);
    // Apply the given label on the graph to on all edges in F (and delete from others)
    void setLabel(const std::string&);

private:
    struct Impl;
    std::unique_ptr<Impl> pimpl_;
};

std::ostream& operator<<(std::ostream& os, const PathFamily& path);
}
