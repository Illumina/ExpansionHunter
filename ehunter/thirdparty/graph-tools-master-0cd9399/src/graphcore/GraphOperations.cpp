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
    Graph reversed(graph.numNodes(), "");

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
