//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "region_spec/RegionBlueprint.hh"

#include <cassert>

using std::string;
using std::vector;

static int getNumberOfRepeats(const string& encoding)
{
    int numBrackets = std::count(encoding.begin(), encoding.end(), '(');
    return numBrackets == 0 ? 1 : numBrackets;
}

RegionBlueprint::RegionBlueprint(
    const string& leftFlank, const string& regionStructureEncoding, const string& rightFlank,
    const vector<string>& repeatIds, const vector<Region>& repeatReferenceRegions,
    const vector<RegionBlueprintComponent::Rarity>& repeatRarities)
{
    assert(repeatIds.size() == std::size_t(getNumberOfRepeats(regionStructureEncoding)));
    assert(repeatReferenceRegions.size() == std::size_t(getNumberOfRepeats(regionStructureEncoding)));
    assert(repeatRarities.size() == std::size_t(getNumberOfRepeats(regionStructureEncoding)));

    components_.emplace_back(
        "LF", leftFlank, RegionBlueprintComponent::Type::kFlank, RegionBlueprintComponent::Rarity::kRare);

    const auto regionStructure = decodeRegionBlueprintSequence(regionStructureEncoding);
    int repeatIndex = 0;

    for (std::size_t index = 0; index != regionStructure.size(); ++index)
    {
        string componentId;
        const auto& sequence = regionStructure[index].sequence;
        const auto type = regionStructure[index].label;
        auto rarity = RegionBlueprintComponent::Rarity::kRare;

        if (type == RegionBlueprintComponent::Type::kRepeat)
        {
            rarity = repeatRarities[repeatIndex];
            componentId = repeatIds[repeatIndex];
        }

        components_.emplace_back(componentId, sequence, type, rarity);

        if (type == RegionBlueprintComponent::Type::kRepeat)
        {
            const auto& referenceRegion = repeatReferenceRegions[repeatIndex];
            components_.back().setReferenceRegion(referenceRegion);
            ++repeatIndex;
        }
    }

    components_.emplace_back(
        "RF", rightFlank, RegionBlueprintComponent::Type::kFlank, RegionBlueprintComponent::Rarity::kRare);
}

std::ostream& operator<<(std::ostream& out, RegionBlueprintComponent::Type componentType)
{
    switch (componentType)
    {
    case RegionBlueprintComponent::Type::kFlank:
        out << "Flank";
        break;
    case RegionBlueprintComponent::Type::kInterruption:
        out << "Interruption";
        break;
    case RegionBlueprintComponent::Type::kRepeat:
        out << "Repeat";
        break;
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, RegionBlueprintComponent::Rarity componentRarity)
{
    switch (componentRarity)
    {
    case RegionBlueprintComponent::Rarity::kRare:
        out << "rare";
        break;
    case RegionBlueprintComponent::Rarity::kCommon:
        out << "common";
        break;
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const RegionBlueprintComponent& component)
{
    out << "[" << component.id() << "," << component.type() << "," << component.rarity() << "," << component.sequence()
        << "]";

    return out;
}

vector<LabeledSequence<RegionBlueprintComponent::Type>> decodeRegionBlueprintSequence(const string& encoding)
{
    vector<LabeledSequence<RegionBlueprintComponent::Type>> sequenceComponents;

    // Encoding uses single-unit format
    if (getNumberOfRepeats(encoding) == 1)
    {
        sequenceComponents.emplace_back(encoding, RegionBlueprintComponent::Type::kRepeat);
    }
    else // Encoding uses multi-unit format
    {
        string currentComponentSequence;
        auto currentComponentType = RegionBlueprintComponent::Type::kRepeat;

        for (char currentSymbol : encoding)
        {
            const bool isCurrentSymbolLeftBracket = currentSymbol == '(';
            const bool isCurrentSymbolRightBracket = currentSymbol == ')';
            const bool isCurrentSymbolBracket = isCurrentSymbolLeftBracket || isCurrentSymbolRightBracket;

            if (isCurrentSymbolBracket && !currentComponentSequence.empty())
            {
                sequenceComponents.emplace_back(currentComponentSequence, currentComponentType);
                currentComponentSequence.clear();
            }

            if (isCurrentSymbolLeftBracket)
            {
                currentComponentType = RegionBlueprintComponent::Type::kRepeat;
            }
            else if (isCurrentSymbolRightBracket)
            {
                currentComponentType = RegionBlueprintComponent::Type::kInterruption;
            }
            else
            {
                currentComponentSequence += currentSymbol;
            }
        }
    }

    return sequenceComponents;
}
