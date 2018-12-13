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

#include "graphcore/PathFamily.hh"

#include <algorithm>
#include <cassert>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include "graphcore/Path.hh"

namespace graphtools
{
struct PathFamily::Impl
{
    explicit Impl(Graph* const graph)
        : graph_(graph)
    {
    }

    Graph* const graph_;
    std::unordered_set<NodeIdPair> edges_;
    std::unordered_set<NodeId> inNodes_;
    std::unordered_set<NodeId> outNodes_;
};

PathFamily::PathFamily(Graph* const graph)
    : pimpl_(new PathFamily::Impl(graph))
{
}

PathFamily::PathFamily(Graph* const graph, const std::string& label)
    : pimpl_(new PathFamily::Impl(graph))
{
    for (const auto& edge : graph->edgesWithLabel(label))
    {
        addEdge(edge.first, edge.second);
    }
}

PathFamily::PathFamily(const PathFamily& other)
    : pimpl_(new Impl(*other.pimpl_))
{
}

PathFamily::PathFamily(PathFamily&& other) noexcept
    : pimpl_(std::move(other.pimpl_))
{
}

PathFamily::~PathFamily() noexcept = default;

const Graph& PathFamily::graph() const { return *(pimpl_->graph_); }
const std::unordered_set<NodeIdPair>& PathFamily::edges() const { return pimpl_->edges_; }

bool PathFamily::operator==(const PathFamily& other) const
{
    return ((edges() == other.edges()) && (pimpl_->graph_ == other.pimpl_->graph_));
}

bool PathFamily::operator!=(const PathFamily& other) const { return !(*this == other); }

bool PathFamily::containsPath(const Path& path) const
{
    int matched(0);
    for (auto start = path.begin(); start != std::prev(path.end()); ++start)
    {
        auto end = std::next(start);
        if (edges().count(NodeIdPair(*start, *end)))
        {
            ++matched;
        }
        else
        {
            if (pimpl_->outNodes_.count(*start) || pimpl_->inNodes_.count(*end))
            {
                return false;
            }
        }
    }
    return matched > 0;
}

bool PathFamily::includes(const PathFamily& other) const
{
    return std::all_of(
        other.edges().begin(), other.edges().end(), [&](const NodeIdPair& e) { return edges().count(e); });
}

void PathFamily::addEdge(NodeId first, NodeId second)
{
    if (!graph().hasEdge(first, second))
    {
        throw std::logic_error("Edges added to path family is not in the graph.");
    }
    pimpl_->edges_.emplace(first, second);
    pimpl_->outNodes_.insert(first);
    pimpl_->inNodes_.insert(second);
}

void PathFamily::setLabel(const std::string& label)
{
    pimpl_->graph_->eraseLabel(label);
    for (const auto& edge : edges())
    {
        pimpl_->graph_->addLabelToEdge(edge.first, edge.second, label);
    }
}

std::ostream& operator<<(std::ostream& os, const PathFamily& path)
{
    os << "{";
    for (const auto& edge : path.edges())
    {
        os << "(" << edge.first << ", " << edge.second << "), ";
    }
    os << "}";
    return os;
}
}
