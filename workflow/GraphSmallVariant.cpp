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

#include "workflow/GraphSmallVariant.hh"

namespace ehunter
{

GraphSmallVariant::GraphSmallVariant(std::shared_ptr<GraphModel> modelPtr, std::vector<graphtools::NodeId> nodeIds)
    : GraphVariant(modelPtr, std::move(nodeIds))
    , alignmentClassifier_(nodeIds_)

{
}

void GraphSmallVariant::process(
    const Read& read, const Alignments& readAligns, const Read& mate, const Alignments& mateAligns)
{
    processRead(read, readAligns);
    processRead(mate, mateAligns);
}

void GraphSmallVariant::processRead(const Read& read, const std::list<graphtools::GraphAlignment>& alignments)
{
    ReadSummaryForSmallVariant smallVariantRead = alignmentClassifier_.classifyRead(read.sequence(), alignments);

    if (smallVariantRead.numAlignments() > 0)
    {
        readSummaries_.push_back(smallVariantRead);
    }

    const auto& smallVariantAlignment = smallVariantRead.alignments().front();

    if (smallVariantAlignment.type() == SmallVariantAlignment::Type::kSpanning)
    {
        if (smallVariantAlignment.nodeId() == SmallVariantAlignmentClassifier::kInvalidNodeId)
        {
            numBypassingReads_ += 1;
        }
        else
        {
            countsOfSpanningReads_.incrementCountOf(smallVariantAlignment.nodeId());
        }
    }

    if (smallVariantAlignment.type() == SmallVariantAlignment::Type::kUpstreamFlanking)
    {
        countsOfReadsFlankingUpstream_.incrementCountOf(smallVariantAlignment.nodeId());
    }

    if (smallVariantAlignment.type() == SmallVariantAlignment::Type::kDownstreamFlanking)
    {
        countsOfReadsFlankingDownstream_.incrementCountOf(smallVariantAlignment.nodeId());
    }
}

int GraphSmallVariant::countReadsSupportingNode(graphtools::NodeId nodeId) const
{
    if (nodeId == SmallVariantAlignmentClassifier::kInvalidNodeId)
    {
        return numBypassingReads_;
    }

    const int numReadsSupportingUpstreamFlank
        = countsOfReadsFlankingUpstream_.countOf(nodeId) + countsOfSpanningReads_.countOf(nodeId);
    const int numReadsSupportingDownstreamFlank
        = countsOfReadsFlankingDownstream_.countOf(nodeId) + countsOfSpanningReads_.countOf(nodeId);

    return (numReadsSupportingUpstreamFlank + numReadsSupportingDownstreamFlank) / 2;
};

}
