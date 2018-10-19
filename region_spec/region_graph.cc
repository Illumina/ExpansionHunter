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

#include "region_spec/region_graph.h"

using graphtools::Graph;
using graphtools::NodeId;
using std::string;
using std::vector;

static int getNumberOfRepeats(const string& encoding)
{
    int numBrackets = std::count(encoding.begin(), encoding.end(), '(');
    return numBrackets == 0 ? 1 : numBrackets;
}

Graph makeRegionGraph(const string& leftFlank, const string& regionStructureEncoding, const string& rightFlank)
{
    const int numRepeats = getNumberOfRepeats(regionStructureEncoding);

    vector<string> repeatIds;
    vector<Region> repeatReferenceRegions;
    vector<RegionBlueprintComponent::Rarity> repeatRarities;

    for (int index = 0; index != numRepeats; ++index)
    {
        repeatIds.push_back("Repeat" + std::to_string(index));
        repeatReferenceRegions.emplace_back("chr", 1, 2);
        repeatRarities.push_back(RegionBlueprintComponent::Rarity::kCommon);
    }

    RegionBlueprint blueprint(
        leftFlank, regionStructureEncoding, rightFlank, repeatIds, repeatReferenceRegions, repeatRarities);

    return makeRegionGraph(blueprint);
}

Graph makeRegionGraph(const RegionBlueprint& blueprint)
{
    Graph graph(blueprint.size());

    const NodeId lastNodeId = blueprint.size() - 1;

    NodeId currentNodeId = 0;
    for (const auto& component : blueprint)
    {
        graph.setNodeSeq(currentNodeId, component.sequence());

        const NodeId previousNodeId = currentNodeId - 1;
        const NodeId nextNodeId = currentNodeId + 1;

        if (currentNodeId != lastNodeId)
        {
            graph.addEdge(currentNodeId, nextNodeId);
        }

        if (component.type() == RegionBlueprintComponent::Type::kRepeat)
        {
            graph.addEdge(currentNodeId, currentNodeId);
            graph.addEdge(previousNodeId, nextNodeId);
        }

        ++currentNodeId;
    }

    return graph;
}
