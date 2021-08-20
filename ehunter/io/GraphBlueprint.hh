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

bool doesFeatureDefineVariant(GraphBlueprintFeatureType featureType);
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
