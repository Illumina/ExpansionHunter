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

#include "graphcore/GraphBuilders.hh"

// needed in gcc 5
#include <cmath>

using std::string;

namespace graphtools
{

Graph makeDeletionGraph(const string& left_flank, const string& deletion, const string& right_flank)
{
    Graph graph(3);

    graph.setNodeSeq(0, left_flank);
    graph.setNodeSeq(1, deletion);
    graph.setNodeSeq(2, right_flank);
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(1, 2);

    return graph;
}

Graph makeSwapGraph(
    const string& left_flank, const string& deletion, const string& insertion, const string& right_flank)
{
    Graph graph(4);

    graph.setNodeSeq(0, left_flank);
    graph.setNodeSeq(1, deletion);
    graph.setNodeSeq(2, insertion);
    graph.setNodeSeq(3, right_flank);
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(1, 3);
    graph.addEdge(2, 3);

    return graph;
}

Graph makeDoubleSwapGraph(
    const string& left_flank, const string& deletion1, const string& insertion1, const string& middle,
    const string& deletion2, const string& insertion2, const string& right_flank)
{
    Graph graph(7);

    graph.setNodeSeq(0, left_flank);
    graph.setNodeSeq(1, deletion1);
    graph.setNodeSeq(2, insertion1);
    graph.setNodeSeq(3, middle);
    graph.setNodeSeq(4, deletion2);
    graph.setNodeSeq(5, insertion2);
    graph.setNodeSeq(6, right_flank);
    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(1, 3);
    graph.addEdge(2, 3);
    graph.addEdge(3, 4);
    graph.addEdge(3, 5);
    graph.addEdge(4, 6);
    graph.addEdge(5, 6);

    return graph;
}

Graph makeLooplessStrGraph(
    size_t read_len, const string& left_flank, const string& repeat_unit, const string& right_flank)
{
    const auto num_repeat_unit_nodes = (size_t)std::ceil(read_len / (double)repeat_unit.length());
    const size_t num_nodes = num_repeat_unit_nodes + 2; // Account for flanks

    Graph graph(num_nodes);

    NodeId right_flank_node_id = static_cast<NodeId>(num_nodes - 1);

    graph.setNodeSeq(0, left_flank);
    graph.setNodeSeq(right_flank_node_id, right_flank);
    graph.addEdge(0, right_flank_node_id);

    for (NodeId node_id = 0; node_id != num_repeat_unit_nodes; ++node_id)
    {
        graph.setNodeSeq(node_id + 1, repeat_unit);
        graph.addEdge(node_id, node_id + 1);
        graph.addEdge(node_id + 1, right_flank_node_id);
    }

    return graph;
}

Graph makeStrGraph(const string& left_flank, const string& repeat_unit, const string& right_flank)
{
    Graph graph(3);
    graph.setNodeSeq(0, left_flank);
    graph.setNodeSeq(1, repeat_unit);
    graph.setNodeSeq(2, right_flank);

    graph.addEdge(0, 1);
    graph.addEdge(0, 2);
    graph.addEdge(1, 1);
    graph.addEdge(1, 2);

    return graph;
}
}
