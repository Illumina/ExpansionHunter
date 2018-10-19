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

    if (prefixQueryLength < prefixReferenceLength)
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

    if (suffixQueryLength < suffixReferenceLength)
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
