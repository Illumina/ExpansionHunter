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

#include "graphcore/Path.hh"
#include "graphcore/PathFamily.hh"
#include "graphcore/PathOperations.hh"

#include <algorithm>
#include <cassert>
#include <list>
#include <map>
#include <set>
#include <vector>

using std::list;
using std::map;
using std::set;
using std::vector;

namespace graphtools
{

std::list<Path> getPathSegmentsForFamily(graphtools::PathFamily const& family)
{
    std::list<Path> segments;
    Graph const& graph = family.graph();

    // sort edges to have their start nodes in topological order
    vector<NodeIdPair> edges{ family.edges().begin(), family.edges().end() };
    std::sort(edges.begin(), edges.end(), [](NodeIdPair const& a, NodeIdPair const& b) -> bool {
        return a.first == b.first ? a.second < b.second : a.first < b.first;
    });

    // compute in and out degree of nodes for subgraph given by path family edges
    map<NodeId, size_t> in_count;
    map<NodeId, size_t> out_count;
    for (const auto& edge : edges)
    {
        out_count[edge.first]++;
        in_count[edge.second]++;
    }

    // concatenate path segments within the family
    for (const auto& edge : edges)
    {
        if (edge.first == edge.second)
        {
            continue;
        }

        bool any_path_extended_by_edge = false;
        for (auto& prefix : segments)
        {
            // only extend when we can do so uniquely
            if (prefix.nodeIds().back() == edge.first && in_count[edge.first] == 1 && out_count[edge.first] == 1)
            {
                prefix.extendEndToIncludeNode(edge.second);
                any_path_extended_by_edge = true;
            }
        }
        if (!any_path_extended_by_edge)
        {
            segments.emplace_back(Path{
                &graph, 0, { edge.first, edge.second }, static_cast<int32_t>(graph.nodeSeq(edge.second).size()) });
        }
    }

    return segments;
}

bool enumeratePathCombinationsInFamily(
    PathFamily const& family, std::list<Path> const& segments, std::list<Path>* paths, size_t maxPaths)
{
    if (paths == nullptr)
    {
        throw std::logic_error("paths cannot be null.");
    }
    paths->clear();

    bool complete = true;
    map<NodeId, set<Path>> segmentsStartingAt;
    map<NodeId, set<Path>> segmentsEndingAt;
    for (const auto& segment : segments)
    {
        segmentsStartingAt[segment.nodeIds().front()].insert(segment);
        segmentsEndingAt[segment.nodeIds().back()].insert(segment);
    }

    bool any_merged = true;
    while (any_merged)
    {
        any_merged = false;
        set<Path> merged_subpaths;
        for (const auto& edge : family.edges())
        {
            auto check_and_merge = [&](std::set<Path> const& prefixes, std::set<Path> const& suffixes) {
                for (const auto& prefix : prefixes)
                {
                    for (const auto& suffix : suffixes)
                    {
                        if (checkPathPrefixSuffixOverlap(prefix, suffix) || checkIfPathsAdjacent(prefix, suffix))
                        {
                            Path segment = mergePaths(prefix, suffix);
                            segmentsStartingAt[segment.nodeIds().front()].insert(segment);
                            segmentsEndingAt[segment.nodeIds().back()].insert(segment);

                            merged_subpaths.insert(prefix);
                            merged_subpaths.insert(suffix);
                            any_merged = true;
                        }
                    }
                }
            };

            check_and_merge(segmentsEndingAt[edge.first], segmentsStartingAt[edge.first]);
            check_and_merge(segmentsEndingAt[edge.second], segmentsStartingAt[edge.second]);
            check_and_merge(segmentsEndingAt[edge.first], segmentsStartingAt[edge.second]);
        }
        for (const auto& path : merged_subpaths)
        {
            segmentsStartingAt[path.nodeIds().front()].erase(path);
            segmentsEndingAt[path.nodeIds().back()].erase(path);
        }

        // check we're not over the maximum count
        size_t count = 0;
        for (const auto& start_list : segmentsStartingAt)
        {
            count += start_list.second.size();
        }
        if (count > maxPaths)
        {
            complete = false;
            break;
        }
    }

    for (const auto& start_list : segmentsStartingAt)
    {
        for (const auto& path : start_list.second)
        {
            paths->push_back(path);
            if (paths->size() > maxPaths)
            {
                complete = false;
                break;
            }
        }
        if (paths->size() > maxPaths)
        {
            break;
        }
    }

    return complete;
}

bool getMaximalPathsForFamily(graphtools::PathFamily const& family, std::list<Path>* paths, size_t maxPaths)
{
    const auto segments = getPathSegmentsForFamily(family);
    return enumeratePathCombinationsInFamily(family, segments, paths, maxPaths);
}

std::map<std::string, graphtools::PathFamily> getPathFamiliesFromGraph(graphtools::Graph& graph)
{
    std::map<std::string, graphtools::PathFamily> families;

    for (const auto& label : graph.allLabels())
    {
        families.insert(std::make_pair(label, PathFamily(&graph, label)));
    }

    return families;
}

graphtools::PathFamily pathToPathFamily(graphtools::Graph& graph, graphtools::Path const& path)
{
    graphtools::PathFamily family(&graph);

    for (size_t i = 1; i < path.numNodes(); ++i)
    {
        family.addEdge(path.getNodeIdByIndex(i - 1), path.getNodeIdByIndex(i));
    }

    return family;
}
}
