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

#include "io/RegionGraph.hh"

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "io/GraphBlueprint.hh"

using graphtools::Graph;
using graphtools::NodeId;
using std::string;
using std::vector;

namespace ehunter
{

int getNumNodes(const GraphBlueprint& blueprint)
{
    int numNodes = 0;
    for (const auto& feature : blueprint)
    {
        numNodes += feature.sequences.size();
    }

    return numNodes;
}

static void setFeatureSequences(const GraphBlueprintFeature& feature, Graph& graph)
{
    assert(feature.nodeIds.size() == feature.sequences.size());
    for (int index = 0; index != static_cast<int>(feature.sequences.size()); ++index)
    {
        graph.setNodeSeq(feature.nodeIds[index], feature.sequences[index]);
    }
}

static void
connectFeatures(const GraphBlueprintFeature& sourceFeature, const GraphBlueprintFeature& sinkFeature, Graph& graph)
{
    for (NodeId nodeIdOfPreviousFeature : sourceFeature.nodeIds)
    {
        for (NodeId nodeIdOfNewFeature : sinkFeature.nodeIds)
        {
            graph.addEdge(nodeIdOfPreviousFeature, nodeIdOfNewFeature);
        }
    }
}

static void setInternalFeatureEdges(const GraphBlueprintFeature& feature, Graph& graph)
{
    const bool isRepeat = feature.type == GraphBlueprintFeatureType::kSkippableRepeat
        || feature.type == GraphBlueprintFeatureType::kUnskippableRepeat;

    if (isRepeat)
    {
        assert(feature.nodeIds.size() == 1);
        const NodeId nodeId = feature.nodeIds.front();
        graph.addEdge(nodeId, nodeId);
    }
}

void setOutgoingFeatureEdges(const GraphBlueprint& blueprint, int index, Graph& graph)
{
    const GraphBlueprintFeature& currentFeature = blueprint[index];
    const GraphBlueprintFeature* downstreamFeaturePtr = &blueprint[++index];

    while (isSkippable(downstreamFeaturePtr->type))
    {
        connectFeatures(currentFeature, *downstreamFeaturePtr, graph);
        downstreamFeaturePtr = &blueprint[++index];
    }

    connectFeatures(currentFeature, *downstreamFeaturePtr, graph);
}

Graph makeRegionGraph(const GraphBlueprint& blueprint, const std::string& locusId)
{
    // Implicit assumptions about the graph structure
    assert(blueprint.front().type == GraphBlueprintFeatureType::kLeftFlank);
    assert(blueprint.back().type == GraphBlueprintFeatureType::kRightFlank);

    Graph graph(getNumNodes(blueprint), locusId);

    for (const auto& feature : blueprint)
    {
        setFeatureSequences(feature, graph);
        setInternalFeatureEdges(feature, graph);
    }

    for (int index = 0; index != static_cast<int>(blueprint.size()) - 1; ++index)
    {
        setOutgoingFeatureEdges(blueprint, index, graph);
    }

    return graph;
}

}
