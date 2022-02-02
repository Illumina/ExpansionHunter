//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
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

#include "graphcore/PathOperations.hh"

#include <cassert>
#include <limits>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

#include "graphutils/BaseMatching.hh"

using std::list;
using std::set;
using std::string;
using std::vector;

namespace graphtools
{

list<Path> extendPathStart(const Path& path, int32_t extension_len)
{
    list<Path> extended_paths;

    const NodeId start_node_id = path.nodeIds().front();

    // Start position gives the maximum extension.
    if (extension_len <= path.startPosition())
    {
        Path extended_path(path);
        extended_path.shiftStartAlongNode(extension_len);
        extended_paths.push_back(extended_path);
    }
    else
    {
        const set<NodeId> pred_node_ids = path.graphRawPtr()->predecessors(start_node_id);
        const int32_t leftover_length = extension_len - path.startPosition();
        for (NodeId pred_node_id : pred_node_ids)
        {
            Path path_with_this_node(path);
            path_with_this_node.extendStartToNode(pred_node_id);
            list<Path> extensions_of_path_with_this_node = extendPathStart(path_with_this_node, leftover_length);

            extended_paths.splice(extended_paths.end(), extensions_of_path_with_this_node);
        }
    }

    return extended_paths;
}

list<Path> extendPathEnd(const Path& path, int32_t extension_len)
{
    list<Path> extended_paths;

    const NodeId end_node_id = path.nodeIds().back();

    const auto end_node_length = static_cast<int32_t>(path.graphRawPtr()->nodeSeq(end_node_id).length());
    const int32_t max_extension_at_end_node = end_node_length - path.endPosition();

    if (extension_len <= max_extension_at_end_node)
    {
        Path extended_path(path);
        extended_path.shiftEndAlongNode(extension_len);
        extended_paths.push_back(extended_path);
    }
    else
    {
        const set<NodeId> succ_node_ids = path.graphRawPtr()->successors(end_node_id);
        const int32_t leftover_length = extension_len - max_extension_at_end_node;
        for (NodeId succ_node_id : succ_node_ids)
        {
            Path path_with_this_node(path);
            path_with_this_node.extendEndToNode(succ_node_id);

            list<Path> extensions_of_path_with_this_node = extendPathEnd(path_with_this_node, leftover_length);
            extended_paths.splice(extended_paths.end(), extensions_of_path_with_this_node);
        }
    }

    return extended_paths;
}

list<Path> extendPath(const Path& path, int32_t start_extension_len, int32_t end_extension_len)
{
    list<Path> extended_paths;
    list<Path> start_extended_paths = extendPathStart(path, start_extension_len);
    for (const Path& start_extended_path : start_extended_paths)
    {
        list<Path> end_extended_paths = extendPathEnd(start_extended_path, end_extension_len);
        extended_paths.splice(extended_paths.end(), end_extended_paths);
    }
    return extended_paths;
}

Path extendPathEndMatching(Path path, const std::string& query, size_t qpos)
{
    const Graph& graph = *path.graphRawPtr();
    size_t pos_in_query = qpos + path.length();
    NodeId node_in_graph = path.nodeIds().back();
    size_t pos_in_node = (size_t)path.endPosition();

    std::vector<NodeId> nodes = path.nodeIds();
    bool moved = true;

    while (moved)
    {
        moved = false;
        const std::string& node_sequence = graph.nodeSeq(node_in_graph);

        while (pos_in_query < query.size() && pos_in_node < node_sequence.size()
               && checkIfReferenceBaseMatchesQueryBase(node_sequence.at(pos_in_node), query.at(pos_in_query)))
        {
            moved = true;
            ++pos_in_node;
            ++pos_in_query;
        }

        if (pos_in_node >= node_sequence.size())
        {
            const auto& successors = graph.successors(node_in_graph);

            size_t num_longest_matches = 0;
            size_t current_longest_match = 0;
            NodeId current_successor = 0;

            size_t successor_min_size = std::numeric_limits<size_t>::max();
            for (auto successor : successors)
            {
                successor_min_size = std::min(successor_min_size, graph.nodeSeq(successor).size());
            }

            for (auto successor : successors)
            {
                size_t pos_in_successor = 0;
                const std::string& successor_sequence = graph.nodeSeq(successor);
                while (pos_in_successor < successor_min_size && pos_in_query + pos_in_successor < query.size()
                       && checkIfReferenceBaseMatchesQueryBase(
                              successor_sequence.at(pos_in_successor), query[pos_in_query + pos_in_successor]))
                {
                    ++pos_in_successor;
                }

                if (pos_in_successor > current_longest_match)
                {
                    current_longest_match = pos_in_successor;
                    current_successor = successor;
                    num_longest_matches = 1;
                }
                else if (pos_in_successor == current_longest_match)
                {
                    ++num_longest_matches;
                }
            }

            if (current_longest_match == 0 || num_longest_matches != 1)
            {
                break;
            }

            nodes.push_back(current_successor);
            pos_in_query += current_longest_match;
            pos_in_node = current_longest_match;
            node_in_graph = current_successor;
            moved = true;
        }
    }
    return Path{ &graph, path.startPosition(), nodes, static_cast<int32_t>(pos_in_node) };
}

Path extendPathStartMatching(Path path, const std::string& query, size_t& pos_in_query)
{
    const Graph& graph = *path.graphRawPtr();
    NodeId node_in_graph = path.nodeIds().front();
    auto pos_in_node = (size_t)path.startPosition();

    std::vector<NodeId> nodes = path.nodeIds();
    bool moved = true;

    while (moved)
    {
        moved = false;

        const std::string& node_sequence = graph.nodeSeq(node_in_graph);

        while (pos_in_query > 0 && pos_in_node > 0
               && checkIfReferenceBaseMatchesQueryBase(node_sequence.at(pos_in_node - 1), query.at(pos_in_query - 1)))
        {
            moved = true;
            --pos_in_node;
            --pos_in_query;
        }

        if (pos_in_node == 0)
        {
            const auto& predecessors = graph.predecessors(node_in_graph);

            size_t num_longest_matches = 0;
            size_t current_longest_match = 0;
            NodeId current_predecessor = 0;

            size_t predecessor_min_size = std::numeric_limits<size_t>::max();
            for (auto predecessor : predecessors)
            {
                predecessor_min_size = std::min(predecessor_min_size, graph.nodeSeq(predecessor).size());
            }
            for (auto predecessor : predecessors)
            {
                const std::string& predecessor_sequence = graph.nodeSeq(predecessor);
                size_t pos_in_predecessor = predecessor_sequence.size();
                size_t match_length = 0;
                while (pos_in_predecessor > (predecessor_sequence.size() - predecessor_min_size)
                       && pos_in_query - match_length > 0
                       && checkIfReferenceBaseMatchesQueryBase(
                              predecessor_sequence.at(pos_in_predecessor - 1), query[pos_in_query - match_length - 1]))
                {
                    --pos_in_predecessor;
                    ++match_length;
                }

                if (match_length > current_longest_match)
                {
                    current_longest_match = match_length;
                    current_predecessor = predecessor;
                    num_longest_matches = 1;
                }
                else if (match_length == current_longest_match)
                {
                    ++num_longest_matches;
                }
            }

            if (current_longest_match == 0
                || num_longest_matches != 1) // also gets the case where current_longest_match == 0
            {
                break;
            }

            nodes.insert(nodes.begin(), current_predecessor);
            pos_in_query -= current_longest_match;
            node_in_graph = current_predecessor;
            pos_in_node = graph.nodeSeq(node_in_graph).size() - current_longest_match;
            moved = true;
        }
    }

    return Path{ &graph, static_cast<int32_t>(pos_in_node), nodes, path.endPosition() };
}

Path extendPathMatching(Path path, const std::string& query, size_t& pos_in_query)
{
    return extendPathStartMatching(extendPathEndMatching(std::move(path), query, pos_in_query), query, pos_in_query);
}

vector<string> splitSequenceByPath(const Path& path, const string& sequence)
{
    if (path.length() != sequence.length())
    {
        throw std::logic_error(
            "Split operation requires that " + path.encode() + " and " + sequence + " have same length");
    }

    vector<string> split_seq;
    size_t cur_position = 0;
    for (int32_t node_index = 0; node_index != (signed)path.numNodes(); ++node_index)
    {
        const size_t length_on_node = path.getNodeOverlapLengthByIndex(static_cast<size_t>(node_index));
        split_seq.push_back(sequence.substr(cur_position, length_on_node));
        cur_position += length_on_node;
    }
    return split_seq;
}

/**
 * Return true if two paths are exactly adjacent
 * (i.e. p1 starts just before p2, or the other way around)
 *
 * @param p1 first path
 * @param p2 second path
 * @return true if the paths are adjacent
 */
bool checkIfPathsAdjacent(Path const& p1, Path const& p2)
{
    // if p1 ends after p2 starts, try the other way around
    if (p1.nodeIds().back() > p2.nodeIds().front())
    {
        return checkIfPathsAdjacent(p2, p1);
    }

    // now p1.nodeIds().back() <= p2.nodeIds().front()
    const auto& graph = *p1.graphRawPtr();

    const auto p1_end_node = p1.nodeIds().back();
    const auto p2_start_node = p2.nodeIds().front();

    if (p1_end_node != p2_start_node && !graph.hasEdge(p1_end_node, p2_start_node))
    {
        return false;
    }

    // now we are left with two cases:
    // p1 ends on node before p2
    // p1 ends on same node as p2

    // adjacent nodes -- check if graph has an edge + if the start / end positions are compatible
    if (p1_end_node != p2_start_node)
    {
        assert(graph.hasEdge(p1_end_node, p2_start_node));
        return p2.startPosition() == 0 && p1.endPosition() == (int32_t)graph.nodeSeq(p1_end_node).size() - 1;
    }

    assert(p1_end_node == p2_start_node);
    return p1.endPosition() + 1 == p2.startPosition();
}

/**
 * Return true if two paths overlap either prefix - suffix or suffix-prefix
 * @param p1 first path
 * @param p2 second path
 * @return true if the paths overlap
 */
bool checkPathPrefixSuffixOverlap(Path const& p1, Path const& p2)
{
    // technically we'd want to check that the two graphs match also
    if (p1.numNodes() == 0 || p2.numNodes() == 0)
    {
        return false;
    }
    if (p1.nodeIds().back() < p2.nodeIds().front() || // p1 ends before p2
        p1.nodeIds().front() > p2.nodeIds().back()) // p1 starts after p2
    {
        return false;
    }

    auto p1_it = p1.nodeIds().begin();
    auto p2_it = p2.nodeIds().begin();

    int shared_nodes = 0;
    while (p1_it != p1.nodeIds().end() && p2_it != p2.nodeIds().end())
    {
        if (*p1_it < *p2_it)
        {
            if (p2_it != p2.nodeIds().begin())
            {
                // paths diverged
                return false;
            }
            // --> ignore non-matching prefix of p1 until paths meet
            ++p1_it;
        }
        else if (*p1_it > *p2_it)
        {
            if (p1_it != p1.nodeIds().begin())
            {
                // paths diverged
                return false;
            }
            // --> ignore non-matching prefix of p2 until paths meet
            ++p2_it;
        }
        else
        { // *p1_it == *p2_it
            // paths have met. They must now match until one of them ends
            ++shared_nodes;
            ++p1_it;
            ++p2_it;
        }
    }

    if (shared_nodes == 0)
    {
        return false;
    }

    // if we only share one node, the paths may not overlap on that node
    if (shared_nodes == 1)
    {
        if (p1_it == p1.nodeIds().end() && p2_it == p2.nodeIds().end())
        {
            if (p1.numNodes() > 1 && p2.numNodes() > 1)
            {
                // if they both have > 1 nodes, they should also share more than one of them;
                // otherwise they cannot both end here
                assert(false);
            }
            else if (p1.numNodes() == 1 && p2.numNodes() > 1)
            {
                // p1 starts here, p2 ends here
                if (p2.endPosition() < p1.startPosition())
                {
                    return false;
                }
            }
            else if (p1.numNodes() > 1 && p2.numNodes() == 1)
            {
                // p2 starts here, p1 ends here
                if (p1.endPosition() < p2.startPosition())
                {
                    return false;
                }
            }
            else if (p1.numNodes() == 1 && p2.numNodes() == 1)
            {
                // both paths on same node, check if they overlap there
                return p1.endPosition() >= p2.startPosition() && p2.endPosition() >= p1.startPosition();
            }
        }
        else if (p1_it != p1.nodeIds().end() && p2_it == p2.nodeIds().end())
        {
            // p2 starts+ends on this node p1 starts -- check that p1 starts before p2 ends
            if (p2.endPosition() < p1.startPosition())
            {
                return false;
            }
        }
        else if (p1_it == p1.nodeIds().end() && p2_it != p2.nodeIds().end())
        {
            // p1 starts+ends on this node p2 starts -- check that p2 starts before p1 ends
            if (p1.endPosition() < p2.startPosition())
            {
                return false;
            }
        }
        else
        {
            // this shouldn't happen. we iterate until one of them reaches end() above
            assert(false);
        }
    }

    return true;
}

/**
 * Paths can be merged if they overlap prefix-suffix / suffix-prefix.
 *
 * @param p1 first path
 * @param p2 second path
 * @return merged path
 */
Path mergePaths(Path const& p1, Path const& p2)
{
    assert(checkPathPrefixSuffixOverlap(p1, p2) || checkIfPathsAdjacent(p1, p2));

    int32_t start = -1;
    int32_t end = -1;
    std::vector<NodeId> nodes;
    auto p1_it = p1.nodeIds().begin();
    auto p2_it = p2.nodeIds().begin();
    while ((p1_it != p1.nodeIds().end()) && (p2_it != p2.nodeIds().end()))
    {
        if (*p1_it < *p2_it)
        {
            if (start < 0)
            {
                start = p1.startPosition();
            }
            nodes.push_back(*p1_it);
            ++p1_it;
        }
        else if (*p1_it > *p2_it)
        {
            if (start < 0)
            {
                start = p2.startPosition();
            }
            nodes.push_back(*p2_it);
            ++p2_it;
        }
        else
        { // *p1_it == *p2_it
            if (start < 0)
            {
                start = std::min(p1.startPosition(), p2.startPosition());
            }
            nodes.push_back(*p1_it);
            ++p1_it;
            ++p2_it;
        }
    }
    if (p1_it == p1.nodeIds().end() && p2_it == p2.nodeIds().end())
    {
        end = std::max(p1.endPosition(), p2.endPosition());
    }
    else if (p1_it != p1.nodeIds().end() && p2_it == p2.nodeIds().end())
    {
        nodes.insert(nodes.end(), p1_it, p1.nodeIds().end());
        end = p1.endPosition();
    }
    else if (p1_it == p1.nodeIds().end() && p2_it != p2.nodeIds().end())
    {
        nodes.insert(nodes.end(), p2_it, p2.nodeIds().end());
        end = p2.endPosition();
    }
    assert(start >= 0 && end >= 0);
    return Path(p1.graphRawPtr(), start, nodes, end);
}

/**
 * Merge a set of paths
 *
 * This will merge paths until none of the resulting paths overlap
 *
 * @param paths a list of paths
 */
void greedyMerge(std::list<Path>& paths)
{
    bool has_merged = true;
    while (has_merged && paths.size() > 1)
    {
        auto path_a = paths.begin();
        has_merged = false;
        while (path_a != paths.end())
        {
            auto path_b = std::next(path_a);
            while (path_b != paths.end())
            {
                if (checkPathPrefixSuffixOverlap(*path_a, *path_b))
                {
                    const Path merged_a_b{ mergePaths(*path_a, *path_b) };
                    paths.erase(path_a);
                    paths.erase(path_b);
                    paths.push_back(merged_a_b);
                    has_merged = true;
                    break;
                }
                ++path_b;
            }
            if (has_merged)
            {
                break;
            }
            ++path_a;
        }
    }
}

/**
 * Merge a set of paths
 *
 * This will merge paths exhaustively, each path is merged with all
 * paths it overlaps until we cannot merge anymore
 *
 * @param paths a list of paths
 */
void exhaustiveMerge(std::list<Path>& paths)
{
    bool has_merged = true;
    while (has_merged && paths.size() > 1)
    {
        auto path_a = paths.begin();
        has_merged = false;

        list<Path> new_paths;
        while (path_a != paths.end())
        {
            auto path_b = paths.begin();
            while (path_b != paths.end())
            {
                if (path_a == path_b)
                {
                    ++path_b;
                    continue;
                }
                if (checkPathPrefixSuffixOverlap(*path_a, *path_b))
                {
                    const Path merged_a_b{ mergePaths(*path_a, *path_b) };
                    const bool a_contained_in_b = merged_a_b.encode() == path_b->encode();
                    const bool b_contained_in_a = merged_a_b.encode() == path_a->encode();
                    const bool a_eq_b = a_contained_in_b && b_contained_in_a;

                    if (a_eq_b)
                    {
                        // keep only one of them
                        new_paths.push_back(*path_b);
                    }
                    else if (a_contained_in_b || b_contained_in_a)
                    {
                        new_paths.push_back(merged_a_b);
                    }
                    else
                    {
                        new_paths.push_back(merged_a_b);
                        new_paths.push_back(*path_a);
                        new_paths.push_back(*path_b);
                    }
                    ++path_b;
                    has_merged = true;
                }
                else
                {
                    new_paths.push_back(*path_b);
                    ++path_b;
                }
            }
            if (has_merged)
            {
                break;
            }
            new_paths.push_back(*path_a);
            ++path_a;
        }
        if (has_merged)
        {
            paths = new_paths;
        }
    }
}

std::list<Path> intersectPaths(Path const& p1, Path const& p2)
{
    std::list<Path> result;

    int32_t start = -1;
    int32_t end = -1;
    std::vector<NodeId> nodes;

    const auto savePath = [&p1, &start, &end, &nodes, &result]() {
        if (!nodes.empty())
        {
            result.emplace_back(p1.graphRawPtr(), start, nodes, end);
            nodes.clear();
            start = -1;
            end = -1;
        }
    };

    auto p1_it = p1.nodeIds().begin();
    auto p2_it = p2.nodeIds().begin();
    while ((p1_it != p1.nodeIds().end()) && (p2_it != p2.nodeIds().end()))
    {
        if (*p1_it < *p2_it)
        {
            savePath();
            ++p1_it;
        }
        else if (*p1_it > *p2_it)
        {
            savePath();
            ++p2_it;
        }
        else
        { // *p1_it == *p2_it
            const auto p1_nodesize = (int32_t)p1.graphRawPtr()->nodeSeq(*p1_it).size();
            const auto p2_nodesize = (int32_t)p2.graphRawPtr()->nodeSeq(*p2_it).size();
            if (p1_nodesize != p2_nodesize)
            {
                throw std::logic_error("Intersecting paths on different graphs is not possible.");
            }

            const int32_t start_p1 = p1_it == p1.nodeIds().begin() ? p1.startPosition() : 0;
            const int32_t start_p2 = p2_it == p2.nodeIds().begin() ? p2.startPosition() : 0;
            const int32_t end_p1 = std::next(p1_it) == p1.nodeIds().end() ? p1.endPosition() : p1_nodesize;
            const int32_t end_p2 = std::next(p2_it) == p2.nodeIds().end() ? p2.endPosition() : p2_nodesize;

            const int32_t start_pos = std::max(start_p1, start_p2);
            const int32_t end_pos = std::min(end_p1, end_p2);

            if (start_pos <= end_pos)
            {
                // start within the node => cannot extend previous matched path
                if (start_pos > 0)
                {
                    savePath();
                }

                if (nodes.empty())
                {
                    // Not sure why cppcheck complains here. start will be read when we call savePath in line 503
                    // cppcheck-suppress unreadVariable
                    start = start_pos;
                }
                else if (!p1.graphRawPtr()->hasEdge(nodes.back(), *p1_it))
                {
                    savePath();
                }

                // Not sure why cppcheck complains here. end will be read when we call savePath in line 503
                // cppcheck-suppress unreadVariable
                end = end_pos;
                nodes.push_back(*p1_it);

                // ends before end of node => cannot combine with match on next node
                if (end_pos + 1 < p1_nodesize)
                {
                    savePath();
                }
            }
            else if (!nodes.empty())
            {
                savePath();
            }

            ++p1_it;
            ++p2_it;
        }
    }
    savePath();

    return result;
}

list<Path> generateSubpathForEachNode(const Path& path)
{
    list<Path> subpaths;

    for (size_t node_index = 0; node_index != path.numNodes(); ++node_index)
    {
        const vector<NodeId> subpath_nodes = { path.getNodeIdByIndex(node_index) };
        const int32_t subpath_start = path.getStartPositionOnNodeByIndex(node_index);
        const int32_t subpath_end = path.getEndPositionOnNodeByIndex(node_index);

        subpaths.emplace_back(path.graphRawPtr(), subpath_start, subpath_nodes, subpath_end);
    }

    return subpaths;
}

bool checkIfBookended(const Path& first_path, const Path& second_path)
{
    const NodeId first_path_end_node = first_path.getNodeIdByIndex(first_path.numNodes() - 1);
    const NodeId second_path_start_node = second_path.getNodeIdByIndex(0);
    const bool are_ends_on_same_node = first_path_end_node == second_path_start_node;
    const bool are_positions_adjacent = first_path.endPosition() == second_path.startPosition();

    if (are_ends_on_same_node && are_positions_adjacent)
    {
        return true;
    }

    const Graph& graph = *first_path.graphRawPtr();
    const auto& successors = graph.successors(first_path_end_node);
    const bool are_ends_on_neighboring_nodes = successors.find(second_path_start_node) != successors.end();
    const auto first_paths_last_node_length = static_cast<int>(graph.nodeSeq(first_path_end_node).length());
    const bool is_first_path_ends_at_node_end = first_path.endPosition() == first_paths_last_node_length;
    const bool is_second_path_starts_at_node_start = second_path.startPosition() == 0;

    if (are_ends_on_neighboring_nodes && is_first_path_ends_at_node_end && is_second_path_starts_at_node_start)
    {
        return true;
    }

    return false;
}

Path concatenatePaths(const Path& first_path, const Path& second_path)
{
    if (!checkIfBookended(first_path, second_path))
    {
        std::ostringstream msg;
        msg << "Paths " << first_path << " and " << second_path << " are not bookended";
        throw std::logic_error(msg.str());
    }

    const NodeId first_path_end_node = first_path.getNodeIdByIndex(first_path.numNodes() - 1);
    const NodeId second_path_start_node = second_path.getNodeIdByIndex(0);
    const bool are_ends_on_same_node = first_path_end_node == second_path_start_node;
    const bool are_positions_adjacent = first_path.endPosition() == second_path.startPosition();

    vector<NodeId> merged_node_ids;
    merged_node_ids.reserve(first_path.nodeIds().size() + second_path.nodeIds().size());
    merged_node_ids.insert(merged_node_ids.end(), first_path.nodeIds().begin(), first_path.nodeIds().end());

    if (are_ends_on_same_node && are_positions_adjacent)
    {
        merged_node_ids.insert(merged_node_ids.end(), second_path.nodeIds().begin() + 1, second_path.nodeIds().end());
    }
    else
    {
        merged_node_ids.insert(merged_node_ids.end(), second_path.nodeIds().begin(), second_path.nodeIds().end());
    }

    return Path(first_path.graphRawPtr(), first_path.startPosition(), merged_node_ids, second_path.endPosition());
}
}
