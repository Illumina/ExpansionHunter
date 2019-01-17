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

#include "input/GraphBlueprint.hh"

#include <algorithm>

using graphtools::NodeId;
using std::string;
using std::vector;

namespace ehunter
{

class TokenizationHelper
{
public:
    TokenizationHelper(const string& regex)
        : regex_(regex)
        , currentSymbolIter_(regex_.begin())
    {
    }

    bool reachedEnd() const { return currentSymbolIter_ == regex_.end(); }
    void advance() { ++currentSymbolIter_; }

    char currentSymbol() const { return *currentSymbolIter_; }
    bool pointingAtBase() const { return kBaseSymbols.find(*currentSymbolIter_) != string::npos; }

    bool pointingAtTokenTerminator() const
    {
        // Last character is always a token terminator
        if (currentSymbolIter_ + 1 == regex_.end())
        {
            return true;
        }

        if (isCountQuantifier(*currentSymbolIter_))
        {
            return true;
        }

        char nextSymbol = *(currentSymbolIter_ + 1);
        if (*currentSymbolIter_ == ')' && !isCountQuantifier(nextSymbol))
        {
            return true;
        }

        if (nextSymbol == '(')
        {
            return true;
        }

        return false;
    }

private:
    bool isCountQuantifier(char symbol) const { return kCountQuantifiers.find(symbol) != string::npos; }

    static const string kBaseSymbols;
    static const string kCountQuantifiers;

    const string& regex_;
    string::const_iterator currentSymbolIter_;
};

const string TokenizationHelper::kBaseSymbols("ACGTBDHKMNSRVWY");
const string TokenizationHelper::kCountQuantifiers("*+?");

vector<string> tokenizeRegex(const string& regex)
{
    vector<string> tokens;
    string token;

    TokenizationHelper tokenizationHelper(regex);
    while (!tokenizationHelper.reachedEnd())
    {
        token += tokenizationHelper.currentSymbol();
        if (tokenizationHelper.pointingAtTokenTerminator())
        {
            tokens.push_back(token);
            token.clear();
        }

        tokenizationHelper.advance();
    }

    return tokens;
}

TokenParser::TokenParser()
    : skippableRepeatRegex_("^\\([ACGTBDHKMNSRVWY]+\\)\\*$")
    , unskippableRepeatRegex_("^\\([ACGTBDHKMNSRVWY]+\\)\\+$")
    , insertionOrDeletionRegex_("^\\([ACGTBDHKMNSRVWY]+\\)\\?$")
    , swapRegex_("^\\([ACGTBDHKMNSRVWY]+\\|[ACGTBDHKMNSRVWY]+\\)$")
    , interruptionRegex_("^[ACGTBDHKMNSRVWY]+$")
{
}

FeatureTypeAndSequences TokenParser::parse(const string& token) const
{
    if (std::regex_match(token, insertionOrDeletionRegex_))
    {
        auto sequence = token.substr(1, token.size() - 3);
        return { GraphBlueprintFeatureType::kInsertionOrDeletion, { sequence } };
    }
    else if (std::regex_match(token, skippableRepeatRegex_))
    {
        auto sequence = token.substr(1, token.size() - 3);
        return { GraphBlueprintFeatureType::kSkippableRepeat, { sequence } };
    }
    else if (std::regex_match(token, unskippableRepeatRegex_))
    {
        auto sequence = token.substr(1, token.size() - 3);
        return { GraphBlueprintFeatureType::kUnskippableRepeat, { sequence } };
    }
    else if (std::regex_match(token, swapRegex_))
    {
        auto sequence = token.substr(1, token.size() - 2);
        int pipePosition = sequence.find('|');
        auto firstAllele = sequence.substr(0, pipePosition);
        auto secondAllele = sequence.substr(pipePosition + 1);
        return { GraphBlueprintFeatureType::kSwap, { firstAllele, secondAllele } };
    }
    else if (std::regex_match(token, interruptionRegex_))
    {
        return { GraphBlueprintFeatureType::kInterruption, { token } };
    }
    else
    {
        throw std::logic_error("Could not parse the token " + token);
    }
}

std::ostream& operator<<(std::ostream& out, GraphBlueprintFeatureType tokenType)
{
    switch (tokenType)
    {
    case GraphBlueprintFeatureType::kLeftFlank:
        out << "GraphBlueprintFeatureType::kLeftFlank";
        break;
    case GraphBlueprintFeatureType::kRightFlank:
        out << "GraphBlueprintFeatureType::kRightFlank";
        break;
    case GraphBlueprintFeatureType::kInsertionOrDeletion:
        out << "GraphBlueprintFeatureType::kInsertionOrDeletion";
        break;
    case GraphBlueprintFeatureType::kInterruption:
        out << "GraphBlueprintFeatureType::kInterruption";
        break;
    case GraphBlueprintFeatureType::kSkippableRepeat:
        out << "GraphBlueprintFeatureType::kSkippableRepeat";
        break;
    case GraphBlueprintFeatureType::kUnskippableRepeat:
        out << "GraphBlueprintFeatureType::kUnskippableRepeat";
        break;
    case GraphBlueprintFeatureType::kSwap:
        out << "GraphBlueprintFeatureType::kSwap";
        break;
    default:
        throw std::logic_error("Encountered unknown token type");
    }

    return out;
}

bool isSkippable(GraphBlueprintFeatureType featureType)
{
    switch (featureType)
    {
    case GraphBlueprintFeatureType::kLeftFlank:
    case GraphBlueprintFeatureType::kRightFlank:
    case GraphBlueprintFeatureType::kInterruption:
    case GraphBlueprintFeatureType::kUnskippableRepeat:
    case GraphBlueprintFeatureType::kSwap:
        return false;
    case GraphBlueprintFeatureType::kSkippableRepeat:
    case GraphBlueprintFeatureType::kInsertionOrDeletion:
        return true;
    default:
        throw std::logic_error("Encountered unrecognized feature blueprint type");
    }
}

GraphBlueprint decodeFeaturesFromRegex(const string& regex)
{
    GraphBlueprint blueprint;

    const auto tokens = tokenizeRegex(regex);
    TokenParser parser;

    NodeId firstUnusedNodeId = 0;
    for (int index = 0; index != static_cast<int>(tokens.size()); ++index)
    {
        const auto& token = tokens[index];
        auto featureTypeAndSequences = parser.parse(token);

        auto featureType = featureTypeAndSequences.first;
        const auto& sequences = featureTypeAndSequences.second;

        if (index == 0)
        {
            if (featureType != GraphBlueprintFeatureType::kInterruption)
            {
                throw std::logic_error("Malformed regular expression " + regex);
            }

            featureType = GraphBlueprintFeatureType::kLeftFlank;
            blueprint.push_back(GraphBlueprintFeature(featureType, sequences, { firstUnusedNodeId }));
            ++firstUnusedNodeId;
        }
        else if (index == static_cast<int>(tokens.size()) - 1)
        {
            if (featureType != GraphBlueprintFeatureType::kInterruption)
            {
                throw std::logic_error("Malformed regular expression " + regex);
            }

            featureType = GraphBlueprintFeatureType::kRightFlank;
            blueprint.push_back(GraphBlueprintFeature(featureType, sequences, { firstUnusedNodeId }));
            ++firstUnusedNodeId;
        }
        else
        {
            vector<NodeId> nodeIds;
            for (int nodeIndex = 0; nodeIndex != static_cast<int>(sequences.size()); ++nodeIndex)
            {
                nodeIds.push_back(firstUnusedNodeId);
                ++firstUnusedNodeId;
            }

            blueprint.push_back(GraphBlueprintFeature(featureType, sequences, nodeIds));
        }
    }

    return blueprint;
}

}
