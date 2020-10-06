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

#include "alignment/OperationsOnAlignments.hh"

#include <cassert>
#include <list>

#include <boost/optional.hpp>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "graphcore/GraphBuilders.hh"

#include "alignment/GreedyAlignmentIntersector.hh"

using graphtools::Alignment;
using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::makeStrGraph;
using graphtools::mergeAlignments;
using graphtools::NodeId;
using graphtools::Path;
using std::list;
using std::string;
using std::to_string;

namespace ehunter
{

GraphAlignment extendWithSoftclip(const GraphAlignment& graphAlignment, int leftSoftclipLen, int rightSoftclipLen)
{
    auto sequenceAlignments = graphAlignment.alignments();

    if (leftSoftclipLen)
    {
        auto& firstAlignment = sequenceAlignments.front();
        int leftSoftclipReferenceStart = firstAlignment.referenceStart();
        Alignment leftSoftclip(leftSoftclipReferenceStart, to_string(leftSoftclipLen) + "S");
        firstAlignment = mergeAlignments(leftSoftclip, firstAlignment);
    }

    if (rightSoftclipLen)
    {
        auto& lastAlignment = sequenceAlignments.back();
        int rightSoftclipReferenceStart = lastAlignment.referenceStart() + lastAlignment.referenceLength();
        Alignment rightSoftclip(rightSoftclipReferenceStart, to_string(rightSoftclipLen) + "S");
        lastAlignment = mergeAlignments(lastAlignment, rightSoftclip);
    }

    return GraphAlignment(graphAlignment.path(), sequenceAlignments);
}

int getNumNonrepeatMatchesUpstream(NodeId nodeId, GraphAlignment alignment)
{
    const list<int> repeatNodeIndexes = alignment.getIndexesOfNode(nodeId);

    if (repeatNodeIndexes.empty())
    {
        return 0;
    }

    const int firstRepeatNodeIndex = repeatNodeIndexes.front();
    int numMatches = 0;

    for (int nodeIndex = 0; nodeIndex != firstRepeatNodeIndex; ++nodeIndex)
    {
        numMatches += alignment[nodeIndex].numMatched();
    }

    return numMatches;
}

int getNumNonrepeatMatchesDownstream(NodeId nodeId, GraphAlignment alignment)
{
    const list<int> repeatNodeIndexes = alignment.getIndexesOfNode(nodeId);

    if (repeatNodeIndexes.empty())
    {
        return 0;
    }

    const int lastRepeatNodeIndex = repeatNodeIndexes.back();
    int numMatches = 0;

    for (int nodeIndex = lastRepeatNodeIndex + 1; nodeIndex != static_cast<int>(alignment.size()); ++nodeIndex)
    {
        numMatches += alignment[nodeIndex].numMatched();
    }

    return numMatches;
}

int scoreAlignmentToNonloopNodes(graphtools::GraphAlignment alignment, LinearAlignmentParameters parameters)
{
    int score = 0;
    const Graph& graph = *alignment.path().graphRawPtr();
    for (int nodeIndex = 0; nodeIndex != (int)alignment.size(); ++nodeIndex)
    {
        NodeId nodeId = alignment.path().getNodeIdByIndex(nodeIndex);
        if (graph.successors(nodeId).find(nodeId) == graph.successors(nodeId).end())
        {
            score += scoreAlignment(
                alignment[nodeIndex], parameters.matchScore, parameters.mismatchScore, parameters.gapOpenScore);
        }
    }

    return score;
}

int countFullOverlaps(NodeId nodeId, GraphAlignment alignment)
{
    const list<int> repeatNodeIndexes = alignment.getIndexesOfNode(nodeId);

    const graphtools::Graph& graph = *alignment.path().graphRawPtr();
    const std::size_t nodeLength = graph.nodeSeq(nodeId).length();

    int numFullOverlaps = 0;
    for (auto nodeIndex : repeatNodeIndexes)
    {
        if (alignment[nodeIndex].referenceLength() == nodeLength)
        {
            ++numFullOverlaps;
        }
    }

    return numFullOverlaps;
}

GraphAlignment computeCanonicalAlignment(const list<GraphAlignment>& alignments)
{
    assert(!alignments.empty());

    if (alignments.size() == 1)
    {
        return alignments.front();
    }

    boost::optional<GraphAlignment> canonicalAlignment = alignments.front();

    for (const auto& alignment : alignments)
    {
        GreedyAlignmentIntersector alignmentIntersector(*canonicalAlignment, alignment);
        canonicalAlignment = alignmentIntersector.intersect();

        if (!canonicalAlignment)
        {
            return alignments.front();
        }
    }

    return *canonicalAlignment;
}

}
