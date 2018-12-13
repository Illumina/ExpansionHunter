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

#include "classification/ClassifierOfAlignmentsToVariant.hh"

#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string/join.hpp>

using graphtools::NodeId;
using std::string;
using std::vector;

namespace ehunter
{

static string encode(const vector<NodeId>& nodeIds)
{
    vector<string> encoding;
    for (auto nodeId : nodeIds)
    {
        encoding.push_back(std::to_string(nodeId));
    }

    return boost::algorithm::join(encoding, ", ");
}

const NodeId ClassifierOfAlignmentsToVariant::kInvalidNodeId = std::numeric_limits<NodeId>::max();

ClassifierOfAlignmentsToVariant::ClassifierOfAlignmentsToVariant(vector<NodeId> targetNodes)
    : targetNodes_(std::move(targetNodes))
{
    if (targetNodes_.empty())
    {
        throw std::logic_error("Cannot create a node bundle without nodes");
    }

    for (int index = 1; index != static_cast<int>(targetNodes_.size()); ++index)
    {
        if (targetNodes_[index] != targetNodes_[index - 1] + 1)
        {
            throw std::logic_error("Bundle " + encode(targetNodes_) + " must contain ordered and consecutive nodes");
        }
    }

    firstBundleNode_ = targetNodes_.front();
    lastBundleNode_ = targetNodes_.back();
}

void ClassifierOfAlignmentsToVariant::classify(const graphtools::GraphAlignment& graphAlignment)
{
    bool pathStartsUpstream = false;
    bool pathEndsDownstream = false;
    bool pathOverlapsTargetNode = false;
    NodeId targetNodeOverlapped = kInvalidNodeId;

    for (auto pathNode : graphAlignment.path().nodeIds())
    {
        if (pathNode < firstBundleNode_)
        {
            pathStartsUpstream = true;
        }
        else if (lastBundleNode_ < pathNode)
        {
            pathEndsDownstream = true;
        }
        else if (firstBundleNode_ <= pathNode && pathNode <= lastBundleNode_)
        {
            pathOverlapsTargetNode = true;
            targetNodeOverlapped = pathNode;
        }
    }

    const bool spanningRead = pathStartsUpstream && pathEndsDownstream;
    const bool upstreamFlankingRead = pathStartsUpstream && pathOverlapsTargetNode;
    const bool downstreamFlankingRead = pathOverlapsTargetNode && pathEndsDownstream;

    if (spanningRead)
    {
        if (targetNodeOverlapped == kInvalidNodeId)
        {
            numBypassingReads_ += 1;
        }
        else
        {
            countsOfSpanningReads_.incrementCountOf(targetNodeOverlapped);
        }
    }
    else if (upstreamFlankingRead)
    {
        countsOfReadsFlankingUpstream_.incrementCountOf(targetNodeOverlapped);
    }
    else if (downstreamFlankingRead)
    {
        countsOfReadsFlankingDownstream_.incrementCountOf(targetNodeOverlapped);
    }
}

}
