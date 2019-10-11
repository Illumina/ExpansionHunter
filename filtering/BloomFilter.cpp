//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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
//

#include "filtering/BloomFilter.hh"

#include <list>

#include "graphcore/Path.hh"
#include "graphcore/PathOperations.hh"
#include "graphutils/SequenceOperations.hh"

using graphtools::expandReferenceSequence;
using graphtools::Graph;
using graphtools::NodeId;
using graphtools::Path;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

BloomFilter::BloomFilter()
    : numBits_(1000000)
    , bits_(numBits_)
{
}

static void addKmerPaths(const list<Path>& kmerPaths, BloomFilter& filter)
{
    for (const Path& kmerPath : kmerPaths)
    {
        if (kmerPath.graphRawPtr()->isSequenceExpansionRequired())
        {
            vector<string> expanded_sequences;
            expandReferenceSequence(kmerPath.seq(), expanded_sequences);
            for (const auto& expansion : expanded_sequences)
            {
                filter.add(expansion);
            }
        }
        else
        {
            filter.add(kmerPath.seq());
        }
    }
}

static void addKmerPathsStartingAtNode(const Graph& graph, NodeId nodeId, int kmerLength, BloomFilter& filter)
{
    for (int pos = 0; pos != static_cast<int>(graph.nodeSeq(nodeId).length()); ++pos)
    {
        Path path(&graph, pos, { nodeId }, pos);
        const auto& paths = extendPath(path, 0, kmerLength);
        addKmerPaths(paths, filter);
    }
}

BloomFilter build(const Graph& graph, int kmerLength)
{
    BloomFilter filter;

    for (NodeId nodeId = 0; nodeId != graph.numNodes(); ++nodeId)
    {
        addKmerPathsStartingAtNode(graph, nodeId, kmerLength, filter);
    }

    return filter;
}

}