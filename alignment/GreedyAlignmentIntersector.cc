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

#include "alignment/GreedyAlignmentIntersector.hh"

#include <algorithm>

using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::NodeId;

boost::optional<GraphAlignment> GreedyAlignmentIntersector::intersect()
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

boost::optional<graphtools::GraphAlignment> GreedyAlignmentIntersector::softclipFirstAlignmentToIntersection() const
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

    return shrankAlignment;
}

bool GreedyAlignmentIntersector::checkIfIntersectionIsConsistent() const
{
    if (nodeIndexOfIntersectionStartOnFirstPath_ == nodeIndexOfIntersectionEndOnFirstPath_)
    {
        return intersectionStart_ < intersectionEnd_;
    }
    return true;
}
