//
// GraphTools library
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
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
