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

#include "alignment/GraphVariantAlignmentStats.hh"

#include <memory>

#include <boost/algorithm/string/join.hpp>

namespace ehunter
{

using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using std::list;
using std::string;
using std::vector;

static string encode(const vector<NodeId>& nodeIds)
{
    vector<string> encoding;
    for (auto nodeId : nodeIds)
    {
        encoding.push_back(std::to_string(nodeId));
    }

    return boost::algorithm::join(encoding, ", ");
}

GraphVariantAlignmentStatsCalculator::GraphVariantAlignmentStatsCalculator(std::vector<graphtools::NodeId> variantNodes)
    : variantNodes_(std::move(variantNodes))
{
    if (variantNodes_.empty())
    {
        throw std::logic_error("Cannot create a node bundle without nodes");
    }

    for (int index = 1; index != static_cast<int>(variantNodes_.size()); ++index)
    {
        if (variantNodes_[index] != variantNodes_[index - 1] + 1)
        {
            throw std::logic_error("Bundle " + encode(variantNodes_) + " must contain ordered and consecutive nodes");
        }
    }

    firstVariantNode_ = variantNodes_.front();
    lastVariantNode_ = variantNodes_.back();
}

void GraphVariantAlignmentStatsCalculator::inspect(const GraphAlignment& alignment)
{
    switch (classify(alignment))
    {
    case Flank::kLeft:
        ++numReadsSpanningLeftBreakpoint_;
        break;
    case Flank::kRight:
        ++numReadsSpanningRightBreakpoint_;
        break;
    case Flank::kBoth:
        ++numReadsSpanningLeftBreakpoint_;
        ++numReadsSpanningRightBreakpoint_;
        break;
    default:
        break;
    }
}

GraphVariantAlignmentStatsCalculator::Flank
GraphVariantAlignmentStatsCalculator::classify(const GraphAlignment& alignment) const
{
    int leftFlankSpan = 0;
    int variantSpan = 0;
    int rightFlankSpan = 0;

    for (int nodeIndex = 0; nodeIndex != static_cast<int>(alignment.path().numNodes()); ++nodeIndex)
    {
        const auto node = alignment.path().getNodeIdByIndex(nodeIndex);
        const auto& alignmentToNode = alignment.alignments().at(nodeIndex);
        int nodeSpan = alignmentToNode.referenceLength();

        if (node < firstVariantNode_)
        {
            leftFlankSpan += nodeSpan;
        }
        else if ((firstVariantNode_ <= node) && (node <= lastVariantNode_))
        {
            variantSpan += nodeSpan;
        }
        else if (lastVariantNode_ < node)
        {
            rightFlankSpan += nodeSpan;
        }
    }

    const bool supportsLeftBreakpoint = (leftFlankSpan >= minSpan_) && (variantSpan + rightFlankSpan >= minSpan_);

    const bool supportsRightBreakpoint = (variantSpan + leftFlankSpan >= minSpan_) && (rightFlankSpan >= minSpan_);

    if (supportsLeftBreakpoint && supportsRightBreakpoint)
    {
        return Flank::kBoth;
    }
    else if (supportsLeftBreakpoint)
    {
        return Flank::kLeft;
    }
    else if (supportsRightBreakpoint)
    {
        return Flank::kRight;
    }
    else
    {
        return Flank::kNeither;
    }
}

GraphVariantAlignmentStats GraphVariantAlignmentStatsCalculator::getStats() const
{
    return { numReadsSpanningLeftBreakpoint_, numReadsSpanningRightBreakpoint_ };
}

std::ostream& operator<<(std::ostream& out, const GraphVariantAlignmentStats& stats)
{
    out << "GraphVariantAlignmentStats(" << stats.numReadsSpanningLeftBreakpoint() << ", "
        << stats.numReadsSpanningRightBreakpoint() << ")";
    return out;
}

}
