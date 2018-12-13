//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Felix Schlesinger <fschlesinger@illumina.com>,
//         Peter Krusche <pkrusche@illumina.com>
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

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
