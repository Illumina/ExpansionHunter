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

#include "region_analysis/SmallVariantAnalyzer.hh"

#include <vector>

using boost::optional;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

void SmallVariantAnalyzer::processMates(
    const Read& read, const list<GraphAlignment>& readAlignments, const Read& mate,
    const std::list<GraphAlignment>& mateAlignments)
{
    processRead(read, readAlignments);
    processRead(mate, mateAlignments);
}

void SmallVariantAnalyzer::processRead(const Read& read, const list<GraphAlignment>& alignments)
{
    ReadSummaryForSmallVariant smallVariantRead = alignmentClassifier_.classifyRead(read.sequence(), alignments);

    if (smallVariantRead.numAlignments() != 1)
    {
        return;
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

int SmallVariantAnalyzer::countReadsSupportingNode(graphtools::NodeId nodeId) const
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
}

std::unique_ptr<VariantFindings> SmallVariantAnalyzer::analyze(const LocusStats& stats) const
{
    NodeId refNode = optionalRefNode_ ? *optionalRefNode_ : SmallVariantAlignmentClassifier::kInvalidNodeId;
    NodeId altNode = SmallVariantAlignmentClassifier::kInvalidNodeId;

    switch (variantSubtype_)
    {
    case VariantSubtype::kInsertion:
        altNode = nodeIds_.front();
        break;
    case VariantSubtype::kDeletion:
        altNode = SmallVariantAlignmentClassifier::kInvalidNodeId;
        break;
    case VariantSubtype::kSwap:
        altNode = (refNode == nodeIds_.front()) ? nodeIds_.back() : nodeIds_.front();
        break;
    case VariantSubtype::kSMN:
        if (refNode != nodeIds_.front())
            throw std::logic_error("Invalid SMN specification");
        altNode = nodeIds_.back();
        break;
    default:
        std::ostringstream encoding;
        encoding << variantSubtype_;
        throw std::logic_error("Invalid small variant subtype: " + encoding.str());
    }

    const int refNodeSupport = countReadsSupportingNode(refNode);
    const int altNodeSupport = countReadsSupportingNode(altNode);

    const double haplotypeDepth = stats.alleleCount() == AlleleCount::kTwo ? stats.depth() / 2 : stats.depth();

    SmallVariantGenotyper smallVariantGenotyper(haplotypeDepth, stats.alleleCount());
    auto genotype = smallVariantGenotyper.genotype(refNodeSupport, altNodeSupport);

    auto refAlleleStatus = allelePresenceChecker_.check(haplotypeDepth, refNodeSupport, altNodeSupport);
    auto altAlleleStatus = allelePresenceChecker_.check(haplotypeDepth, altNodeSupport, refNodeSupport);

    return std::unique_ptr<VariantFindings>(
        new SmallVariantFindings(refNodeSupport, altNodeSupport, refAlleleStatus, altAlleleStatus, genotype));
}

}
