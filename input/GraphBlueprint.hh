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

#pragma once

#include <iostream>
#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include "graphcore/Graph.hh"

namespace ehunter
{

enum class GraphBlueprintFeatureType
{
    kLeftFlank,
    kRightFlank,
    kSkippableRepeat,
    kUnskippableRepeat,
    kInsertionOrDeletion,
    kSwap,
    kInterruption
};

bool isSkippable(GraphBlueprintFeatureType featureType);

using FeatureTypeAndSequences = std::pair<GraphBlueprintFeatureType, std::vector<std::string>>;

std::vector<std::string> tokenizeRegex(const std::string& regex);

class TokenParser
{
public:
    TokenParser();

    FeatureTypeAndSequences parse(const std::string& token) const;

private:
    std::regex skippableRepeatRegex_;
    std::regex unskippableRepeatRegex_;
    std::regex insertionOrDeletionRegex_;
    std::regex swapRegex_;
    std::regex interruptionRegex_;
};

struct GraphBlueprintFeature
{
    GraphBlueprintFeature(
        GraphBlueprintFeatureType type, std::vector<std::string> sequences, std::vector<graphtools::NodeId> nodeIds)
        : type(std::move(type))
        , sequences(std::move(sequences))
        , nodeIds(std::move(nodeIds))
    {
    }
    GraphBlueprintFeatureType type;
    std::vector<std::string> sequences;
    std::vector<graphtools::NodeId> nodeIds;
};

using GraphBlueprint = std::vector<GraphBlueprintFeature>;

GraphBlueprint decodeFeaturesFromRegex(const std::string& regex);

std::ostream& operator<<(std::ostream& out, GraphBlueprintFeatureType tokenType);

}
