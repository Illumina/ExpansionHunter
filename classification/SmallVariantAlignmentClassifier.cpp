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

#include "classification/SmallVariantAlignmentClassifier.hh"

#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string/join.hpp>

using boost::optional;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using std::list;
using std::set;
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

const NodeId SmallVariantAlignmentClassifier::kInvalidNodeId = std::numeric_limits<NodeId>::max();

SmallVariantAlignmentClassifier::SmallVariantAlignmentClassifier(vector<NodeId> targetNodes)
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

ReadSummaryForSmallVariant SmallVariantAlignmentClassifier::classifyRead(
    const std::string& read, const std::list<graphtools::GraphAlignment>& alignments)
{
    ReadSummaryForSmallVariant smallVariantRead(read.length());
    for (const auto& alignment : alignments)
    {
        optional<SmallVariantAlignment> optionalSummary = classify(alignment);
        if (optionalSummary)
        {
            smallVariantRead.addAlignment(*optionalSummary);
        }
    }

    return smallVariantRead;
}

optional<SmallVariantAlignment> SmallVariantAlignmentClassifier::classify(const GraphAlignment& alignment)
{
    bool pathStartsUpstream = false;
    bool pathEndsDownstream = false;
    bool pathOverlapsTargetNode = false;
    NodeId targetNodeOverlapped = kInvalidNodeId;

    for (auto pathNode : alignment.path().nodeIds())
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

    const int score = scoreAlignment(alignment);

    if (spanningRead)
    {
        return SmallVariantAlignment(targetNodeOverlapped, SmallVariantAlignment::Type::kSpanning, score);
    }
    else if (upstreamFlankingRead)
    {
        return SmallVariantAlignment(targetNodeOverlapped, SmallVariantAlignment::Type::kUpstreamFlanking, score);
    }
    else if (downstreamFlankingRead)
    {
        return SmallVariantAlignment(targetNodeOverlapped, SmallVariantAlignment::Type::kDownstreamFlanking, score);
    }
    else
    {
        return boost::none;
    }
}
}
