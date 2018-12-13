//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "input/RegionGraph.hh"

#include <algorithm>
#include <cassert>
#include <utility>

#include "input/GraphBlueprint.hh"

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

Graph makeRegionGraph(const GraphBlueprint& blueprint)
{
    // Implicit assumptions about the graph structure
    assert(blueprint.front().type == GraphBlueprintFeatureType::kLeftFlank);
    assert(blueprint.back().type == GraphBlueprintFeatureType::kRightFlank);

    Graph graph(getNumNodes(blueprint));

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
