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

#include "alignment/GreedyAlignmentIntersector.hh"

#include <algorithm>

#include "graphalign/GraphAlignmentOperations.hh"

namespace ehunter
{

using boost::optional;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::isLocalAlignment;
using graphtools::NodeId;

optional<GraphAlignment> GreedyAlignmentIntersector::intersect()
{
    initialize();

    if (!tryAdvancingIndexesToCommonNode())
    {
        return boost::optional<GraphAlignment>();
    }

    if (checkIfCommonNodeIsLoop())
    {
        advanceIndexesToMatchRemainingIterations();
    }

    advanceIndexesToLastCommonNode();
    computeIntersectionEnds();

    if (!checkIfIntersectionIsConsistent())
    {
        return boost::optional<GraphAlignment>();
    }

    return softclipFirstAlignmentToIntersection();
}

void GreedyAlignmentIntersector::initialize()
{
    nodeIndexOfIntersectionStartOnFirstPath_ = 0;
    nodeIndexOfIntersectionStartOnSecondPath_ = 0;
}

bool GreedyAlignmentIntersector::checkIfAlignmentEndReached(int firstPathIndex, int secondPathIndex)
{
    return firstPathIndex == (int)firstAlignment_.size() || secondPathIndex == (int)secondAlignment_.size();
}

bool GreedyAlignmentIntersector::tryAdvancingIndexesToCommonNode()
{
    while (!checkIfAlignmentEndReached(
        nodeIndexOfIntersectionStartOnFirstPath_, nodeIndexOfIntersectionStartOnSecondPath_))
    {
        const NodeId firstPathNode = firstPath_.getNodeIdByIndex(nodeIndexOfIntersectionStartOnFirstPath_);
        const NodeId secondPathNode = secondPath_.getNodeIdByIndex(nodeIndexOfIntersectionStartOnSecondPath_);

        if (firstPathNode < secondPathNode)
        {
            ++nodeIndexOfIntersectionStartOnFirstPath_;
        }
        else if (secondPathNode < firstPathNode)
        {
            ++nodeIndexOfIntersectionStartOnSecondPath_;
        }
        else
        {
            break;
        }
    }

    return !checkIfAlignmentEndReached(
        nodeIndexOfIntersectionStartOnFirstPath_, nodeIndexOfIntersectionStartOnSecondPath_);
}

bool GreedyAlignmentIntersector::checkIfCommonNodeIsLoop() const
{
    const NodeId firstPathNode = firstPath_.getNodeIdByIndex(nodeIndexOfIntersectionStartOnFirstPath_);
    assert(firstPathNode == secondPath_.getNodeIdByIndex(nodeIndexOfIntersectionStartOnSecondPath_));

    const Graph& graph = *firstPath_.graphRawPtr();
    return graph.hasEdge(firstPathNode, firstPathNode);
}

void GreedyAlignmentIntersector::advanceIndexesToMatchRemainingIterations()
{
    const NodeId loopNodeId = firstPath_.getNodeIdByIndex(nodeIndexOfIntersectionStartOnFirstPath_);
    // TODO: Move getIndexesOfNode to path class;
    const int numIterationsMadeByFirstPath = firstAlignment_.getIndexesOfNode(loopNodeId).size();
    const int numIterationsMadeBySecondPath = secondAlignment_.getIndexesOfNode(loopNodeId).size();

    if (numIterationsMadeByFirstPath < numIterationsMadeBySecondPath)
    {
        nodeIndexOfIntersectionStartOnSecondPath_ += numIterationsMadeBySecondPath - numIterationsMadeByFirstPath;
    }
    else if (numIterationsMadeBySecondPath < numIterationsMadeByFirstPath)
    {
        nodeIndexOfIntersectionStartOnFirstPath_ += numIterationsMadeByFirstPath - numIterationsMadeBySecondPath;
    }
}

void GreedyAlignmentIntersector::advanceIndexesToLastCommonNode()
{
    nodeIndexOfIntersectionEndOnFirstPath_ = nodeIndexOfIntersectionStartOnFirstPath_;
    nodeIndexOfIntersectionEndOnSecondPath_ = nodeIndexOfIntersectionStartOnSecondPath_;

    while (!checkIfAlignmentEndReached(
        nodeIndexOfIntersectionEndOnFirstPath_ + 1, nodeIndexOfIntersectionEndOnSecondPath_ + 1))
    {
        const NodeId firstPathNode = firstPath_.getNodeIdByIndex(nodeIndexOfIntersectionEndOnFirstPath_ + 1);
        const NodeId secondPathNode = secondPath_.getNodeIdByIndex(nodeIndexOfIntersectionEndOnSecondPath_ + 1);

        if (firstPathNode == secondPathNode)
        {
            ++nodeIndexOfIntersectionEndOnFirstPath_;
            ++nodeIndexOfIntersectionEndOnSecondPath_;
        }
        else
        {
            break;
        }
    }
}

void GreedyAlignmentIntersector::computeIntersectionEnds()
{
    const int firstPathStartPosition
        = firstPath_.getStartPositionOnNodeByIndex(nodeIndexOfIntersectionStartOnFirstPath_);
    const int secondPathStartPosition
        = secondPath_.getStartPositionOnNodeByIndex(nodeIndexOfIntersectionStartOnSecondPath_);

    intersectionStart_ = std::max(firstPathStartPosition, secondPathStartPosition);

    const int firstPathEndPosition = firstPath_.getEndPositionOnNodeByIndex(nodeIndexOfIntersectionEndOnFirstPath_);
    const int secondPathEndPosition = secondPath_.getEndPositionOnNodeByIndex(nodeIndexOfIntersectionEndOnSecondPath_);

    intersectionEnd_ = std::min(firstPathEndPosition, secondPathEndPosition);
}

optional<GraphAlignment> GreedyAlignmentIntersector::softclipFirstAlignmentToIntersection() const
{
    GraphAlignment shrankAlignment = firstAlignment_;

    int leftoverPrefixReferenceLength = 0;
    for (int nodeIndex = 0; nodeIndex != nodeIndexOfIntersectionStartOnFirstPath_; ++nodeIndex)
    {
        leftoverPrefixReferenceLength += firstPath_.getNodeOverlapLengthByIndex(nodeIndex);
    }

    const int originalStartPosition
        = firstPath_.getStartPositionOnNodeByIndex(nodeIndexOfIntersectionStartOnFirstPath_);
    leftoverPrefixReferenceLength += intersectionStart_ - originalStartPosition;

    if (leftoverPrefixReferenceLength)
    {
        shrankAlignment.shrinkStart(leftoverPrefixReferenceLength);
    }

    const int numNodes = firstPath_.numNodes();
    int leftoverSuffixReferenceLength = 0;
    for (int nodeIndex = nodeIndexOfIntersectionEndOnFirstPath_ + 1; nodeIndex != numNodes; ++nodeIndex)
    {
        leftoverSuffixReferenceLength += firstPath_.getNodeOverlapLengthByIndex(nodeIndex);
    }

    const int originalEndPosition = firstPath_.getEndPositionOnNodeByIndex(nodeIndexOfIntersectionEndOnFirstPath_);
    leftoverSuffixReferenceLength += originalEndPosition - intersectionEnd_;

    if (leftoverSuffixReferenceLength)
    {
        shrankAlignment.shrinkEnd(leftoverSuffixReferenceLength);
    }

    // Note that shrankAlignment may not always be a local alignment (that is an alignment that starts and ends with a
    // match possibly flanked by soft clips). Hence an explicit check below is required.
    return isLocalAlignment(shrankAlignment) ? shrankAlignment : optional<GraphAlignment>();
}

bool GreedyAlignmentIntersector::checkIfIntersectionIsConsistent() const
{
    if (nodeIndexOfIntersectionStartOnFirstPath_ == nodeIndexOfIntersectionEndOnFirstPath_)
    {
        return intersectionStart_ < intersectionEnd_;
    }
    return true;
}

}
