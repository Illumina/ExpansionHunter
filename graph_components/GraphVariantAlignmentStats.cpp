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

#include "graph_components/GraphVariantAlignmentStats.hh"

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

void GraphVariantAlignmentStatsCalculator::inspect(const list<GraphAlignment>& alignments)
{
    for (const auto& alignment : alignments)
    {
        switch (classify(alignment))
        {
        case Flank::kLeft:
            ++numReadsOverlappingLeftBreakpoint_;
            return;
        case Flank::kRight:
            ++numReadsOverlappingRightBreakpoint_;
            return;
        case Flank::kBoth:
            ++numReadsOverlappingLeftBreakpoint_;
            ++numReadsOverlappingRightBreakpoint_;
            return;
        default:
            break;
        }
    }
}

GraphVariantAlignmentStatsCalculator::Flank
GraphVariantAlignmentStatsCalculator::classify(const GraphAlignment& alignment) const
{
    int numLeftFlankMatches = 0;
    int numStrMatches = 0;
    int numRightFlankMatches = 0;

    for (int nodeIndex = 0; nodeIndex != static_cast<int>(alignment.path().numNodes()); ++nodeIndex)
    {
        const auto node = alignment.path().getNodeIdByIndex(nodeIndex);
        const auto& alignmentToNode = alignment.alignments().at(nodeIndex);
        int numMatches = alignmentToNode.numMatched();

        if (node < firstVariantNode_)
        {
            numLeftFlankMatches += numMatches;
        }
        else if ((firstVariantNode_ <= node) && (node <= lastVariantNode_))
        {
            numStrMatches += numMatches;
        }
        else if (lastVariantNode_ < node)
        {
            numRightFlankMatches += numMatches;
        }
    }

    const bool supportsLeftBreakpoint
        = (numLeftFlankMatches >= minMatch_) && (numStrMatches + numRightFlankMatches >= minMatch_);

    const bool supportsRightBreakpoint
        = (numStrMatches + numLeftFlankMatches >= minMatch_) && (numRightFlankMatches >= minMatch_);

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

GraphVariantAlignmentStats GraphVariantAlignmentStatsCalculator::getStats(int readLength) const
{
    double leftBreakpointCoverage = computeBreakpointCoverage(numReadsOverlappingLeftBreakpoint_, readLength);
    double rightBreakpointCoverage = computeBreakpointCoverage(numReadsOverlappingRightBreakpoint_, readLength);
    return { leftBreakpointCoverage, rightBreakpointCoverage };
}

double GraphVariantAlignmentStatsCalculator::computeBreakpointCoverage(int numReads, int readLength) const
{
    return static_cast<double>(numReads * readLength) / (readLength - 2 * minMatch_);
}

std::ostream& operator<<(std::ostream& out, const GraphVariantAlignmentStats& stats)
{
    out << "StrAlignmentStats(" << stats.leftBreakpointCoverage() << ", " << stats.rightBreakpointCoverage() << ")";
    return out;
}

}