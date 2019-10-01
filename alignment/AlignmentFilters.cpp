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

#include "alignment/AlignmentFilters.hh"

#include <algorithm>
#include <list>
#include <vector>

#include "graphalign/GaplessAligner.hh"
#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphcore/PathOperations.hh"

#include "alignment/OperationsOnAlignments.hh"

using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::Operation;
using graphtools::OperationType;
using graphtools::Path;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

bool checkIfComesFromGraphLocus(
    const list<GraphAlignment>& readAlignments, const list<GraphAlignment>& mateAlignments,
    int kMinNonRepeatAlignmentScore)
{
    int maxReadScore = 0;
    for (const auto& alignment : readAlignments)
    {
        maxReadScore = std::max(maxReadScore, scoreAlignmentToNonloopNodes(alignment));
    }

    int maxMateScore = 0;
    for (const auto& alignment : mateAlignments)
    {
        maxMateScore = std::max(maxMateScore, scoreAlignmentToNonloopNodes(alignment));
    }

    return maxReadScore + maxMateScore >= kMinNonRepeatAlignmentScore;
}

bool checkIfUpstreamAlignmentIsGood(NodeId nodeId, const GraphAlignment& alignment)
{
    const list<int> repeatNodeIndexes = alignment.getIndexesOfNode(nodeId);

    if (repeatNodeIndexes.empty())
    {
        return false;
    }

    const int firstRepeatNodeIndex = repeatNodeIndexes.front();
    int score = 0;
    LinearAlignmentParameters parameters;
    for (int nodeIndex = 0; nodeIndex != firstRepeatNodeIndex; ++nodeIndex)
    {
        score += scoreAlignment(
            alignment[nodeIndex], parameters.matchScore, parameters.mismatchScore, parameters.gapOpenScore);
    }

    const int kScoreCutoff = parameters.matchScore * 8;

    return score >= kScoreCutoff;
}

bool checkIfDownstreamAlignmentIsGood(NodeId nodeId, const GraphAlignment& alignment)
{
    const list<int> repeatNodeIndexes = alignment.getIndexesOfNode(nodeId);

    if (repeatNodeIndexes.empty())
    {
        return false;
    }

    const int lastRepeatNodeIndex = repeatNodeIndexes.back();
    int score = 0;
    LinearAlignmentParameters parameters;
    for (int nodeIndex = lastRepeatNodeIndex + 1; nodeIndex != static_cast<int>(alignment.size()); ++nodeIndex)
    {
        score += scoreAlignment(
            alignment[nodeIndex], parameters.matchScore, parameters.mismatchScore, parameters.gapOpenScore);
    }

    const int kScoreCutoff = parameters.matchScore * 8;

    return score >= kScoreCutoff;
}

bool checkIfPassesAlignmentFilters(const GraphAlignment& alignment)
{
    const Operation& firstOperation = alignment.alignments().front().operations().front();
    const int frontSoftclipLen = firstOperation.type() == OperationType::kSoftclip ? firstOperation.queryLength() : 0;

    const Operation& lastOperation = alignment.alignments().back().operations().back();
    const int backSoftclipLen = lastOperation.type() == OperationType::kSoftclip ? lastOperation.queryLength() : 0;

    const int clippedQueryLength = alignment.queryLength() - frontSoftclipLen - backSoftclipLen;
    const int referenceLength = alignment.referenceLength();

    const int percentQueryMatches = (100 * alignment.numMatches()) / clippedQueryLength;
    const int percentReferenceMatches = (100 * alignment.numMatches()) / referenceLength;

    return percentQueryMatches >= 80 && percentReferenceMatches >= 80;
}

}
