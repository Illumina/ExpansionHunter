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
