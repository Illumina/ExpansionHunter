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

#include "alignment/AlignmentTweakers.hh"

#include <list>

#include "graphalign/GaplessAligner.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphalign/LinearAlignmentParameters.hh"
#include "graphcore/PathOperations.hh"

using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::Path;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

static void shrinkPrefixUntilNodeBoundary(Path& path, int maxShrinkLength)
{
    int accumulatedLength = 0;

    while (path.numNodes() > 1 && accumulatedLength + (int)path.getNodeOverlapLengthByIndex(0) <= maxShrinkLength)
    {
        accumulatedLength += path.getNodeOverlapLengthByIndex(0);
        path.removeStartNode();
    }
}

static void shrinkSuffixUntilNodeBoundary(Path& path, int maxShrinkLength)
{
    int accumulatedLength = 0;
    int nodeIndex = path.numNodes() - 1;

    while (path.numNodes() > 1
           && accumulatedLength + (int)path.getNodeOverlapLengthByIndex(nodeIndex) <= maxShrinkLength)
    {
        accumulatedLength += path.getNodeOverlapLengthByIndex(nodeIndex);
        path.removeEndNode();
        --nodeIndex;
    }
}

static list<Path> computeAlternatePrefixes(Path path, int endLength)
{
    path.shrinkEndBy(path.length());
    return extendPathStart(path, endLength);
}

static list<Path> computeAlternateSuffixes(Path path, int endLength)
{
    path.shrinkStartBy(path.length());
    return extendPathEnd(path, endLength);
}

static list<Path> getHighScoringPaths(const list<Path>& paths, const string& query, int lowScoreCutoff)
{
    LinearAlignmentParameters parameters;
    list<Path> highScoringPaths;

    for (const auto& path : paths)
    {
        const auto alignment = graphtools::alignWithoutGaps(0, path.seq(), query);
        const int score
            = scoreAlignment(alignment, parameters.matchScore, parameters.mismatchScore, parameters.gapOpenScore);

        if (score >= lowScoreCutoff)
        {
            highScoringPaths.push_back(path);
        }
    }

    return highScoringPaths;
}

static std::size_t getSmallestNumberOfNodes(const list<Path>& paths)
{
    assert(!paths.empty());

    std::size_t minNodes = paths.front().numNodes();
    for (const auto& path : paths)
    {
        minNodes = std::min(minNodes, path.numNodes());
    }
    return minNodes;
}

static int computeCommonPrefixLength(const list<Path>& paths)
{
    const int numNodes = getSmallestNumberOfNodes(paths);
    const Path& firstPath = paths.front();

    int prefixLength = 0;

    for (int nodeIndex = 0; nodeIndex != numNodes; ++nodeIndex)
    {
        const int firstPathNodeId = firstPath.getNodeIdByIndex(nodeIndex);

        for (const Path& path : paths)
        {
            int pathNodeId = path.getNodeIdByIndex(nodeIndex);
            if (firstPathNodeId != pathNodeId)
            {
                return prefixLength;
            }
        }

        prefixLength += firstPath.getNodeOverlapLengthByIndex(nodeIndex);
    }

    return prefixLength;
}

static int computeCommonSuffixLength(const list<Path>& paths)
{
    const int numNodes = getSmallestNumberOfNodes(paths);
    const Path& firstPath = paths.front();

    int suffixLength = 0;

    for (int nodeIndex = 0; nodeIndex != numNodes; ++nodeIndex)
    {
        const int firstPathReverseNodeIndex = firstPath.numNodes() - nodeIndex - 1;
        const int firstPathNodeId = firstPath.getNodeIdByIndex(firstPathReverseNodeIndex);

        for (const Path& path : paths)
        {
            const int reverseNodeIndex = path.numNodes() - nodeIndex - 1;
            int pathNodeId = path.getNodeIdByIndex(reverseNodeIndex);
            if (firstPathNodeId != pathNodeId)
            {
                return suffixLength;
            }
        }

        suffixLength += firstPath.getNodeOverlapLengthByIndex(firstPathReverseNodeIndex);
    }

    return suffixLength;
}

static int computeQueryLengthUpToNode(const GraphAlignment& alignment, int terminalNodeIndex)
{
    assert(terminalNodeIndex <= (int)alignment.size());

    int queryLength = 0;
    for (int nodeIndex = 0; nodeIndex != terminalNodeIndex; ++nodeIndex)
    {
        queryLength += alignment[nodeIndex].queryLength();
    }

    return queryLength;
}

void shrinkUncertainPrefix(int referenceLength, const string& query, GraphAlignment& alignment)
{
    Path shrunkPath = alignment.path();
    shrinkPrefixUntilNodeBoundary(shrunkPath, referenceLength);
    const int prefixReferenceLength = alignment.referenceLength() - shrunkPath.length();

    if (prefixReferenceLength == 0)
    {
        return;
    }

    const int numPrefixNodesRemoved = alignment.path().numNodes() - shrunkPath.numNodes();
    const int prefixQueryLength = computeQueryLengthUpToNode(alignment, numPrefixNodesRemoved);

    if (prefixQueryLength != prefixReferenceLength)
    {
        alignment.shrinkStart(prefixReferenceLength);
        return;
    }

    string trimmedQueryPrefix = query.substr(0, prefixQueryLength);
    trimmedQueryPrefix
        = trimmedQueryPrefix.substr(trimmedQueryPrefix.length() - prefixReferenceLength, prefixReferenceLength);

    const list<Path> alternatePrefixes = computeAlternatePrefixes(shrunkPath, prefixReferenceLength);
    assert(!alternatePrefixes.empty());

    LinearAlignmentParameters parameters;
    const int lowScoreCutoff_ = (prefixReferenceLength / 2) * parameters.matchScore;
    const list<Path> highScoringPrefixes = getHighScoringPaths(alternatePrefixes, trimmedQueryPrefix, lowScoreCutoff_);

    if (highScoringPrefixes.empty())
    {
        alignment.shrinkStart(prefixReferenceLength);
        return;
    }

    const int lengthSharedByPrefixes = computeCommonSuffixLength(highScoringPrefixes);

    alignment.shrinkStart(prefixReferenceLength - lengthSharedByPrefixes);
}

void shrinkUncertainSuffix(int referenceLength, const string& query, GraphAlignment& alignment)
{
    Path shrunkPath = alignment.path();
    shrinkSuffixUntilNodeBoundary(shrunkPath, referenceLength);
    const int suffixReferenceLength = alignment.referenceLength() - shrunkPath.length();

    if (suffixReferenceLength == 0)
    {
        return;
    }

    const int prefixQueryLength = computeQueryLengthUpToNode(alignment, shrunkPath.numNodes());
    const int suffixQueryLength = alignment.queryLength() - prefixQueryLength;

    if (suffixQueryLength != suffixReferenceLength)
    {
        alignment.shrinkEnd(suffixReferenceLength);
        return;
    }

    const string trimmedQuerySuffix = query.substr(prefixQueryLength, suffixReferenceLength);

    const list<Path> alternateSuffixes = computeAlternateSuffixes(shrunkPath, suffixReferenceLength);
    assert(!alternateSuffixes.empty());

    LinearAlignmentParameters parameters;
    const int lowScoreCutoff_ = (suffixReferenceLength / 2) * parameters.matchScore;
    const list<Path> highScoringSuffixes = getHighScoringPaths(alternateSuffixes, trimmedQuerySuffix, lowScoreCutoff_);

    if (highScoringSuffixes.empty())
    {
        alignment.shrinkEnd(suffixReferenceLength);
        return;
    }

    const int lengthSharedBySuffixes = computeCommonPrefixLength(highScoringSuffixes);

    alignment.shrinkEnd(suffixReferenceLength - lengthSharedBySuffixes);
}

}
