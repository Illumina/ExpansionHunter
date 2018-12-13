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

#include "graphcore/GraphOperations.hh"

#include "graphutils/SequenceOperations.hh"

namespace graphtools
{
/**
 * Reverse (and optionally complement) a graph
 * @param graph the graph to reverse
 * @param complement true to also reverse-complement all sequences
 */
Graph reverseGraph(Graph const& graph, bool complement)
{
    Graph reversed(graph.numNodes(), "", graph.isSequenceExpansionRequired());

    for (NodeId node_id = 0; node_id != graph.numNodes(); ++node_id)
    {
        reversed.setNodeSeq(
            static_cast<NodeId>(graph.numNodes() - 1 - node_id),
            complement ? reverseComplement(graph.nodeSeq(node_id)) : reverseString(graph.nodeSeq(node_id)));

        for (auto succ : graph.successors(node_id))
        {
            const auto to_new = static_cast<NodeId>(graph.numNodes() - 1 - node_id);
            const auto from_new = static_cast<NodeId>(graph.numNodes() - 1 - succ);
            reversed.addEdge(from_new, to_new);
            for (const auto& label : graph.edgeLabels(node_id, succ))
            {
                reversed.addLabelToEdge(from_new, to_new, label);
            }
        }
    }
    return reversed;
}
}
