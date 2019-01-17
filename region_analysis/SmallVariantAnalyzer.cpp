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

#include "region_analysis/SmallVariantAnalyzer.hh"

#include <vector>

#include <boost/optional.hpp>

using boost::optional;
using graphtools::NodeId;
using std::string;
using std::vector;

namespace ehunter
{

void SmallVariantAnalyzer::processMates(
    const reads::Read& /*read*/, const graphtools::GraphAlignment& readAlignment, const reads::Read& /*mate*/,
    const graphtools::GraphAlignment& mateAlignment)
{
    alignmentClassifier_.classify(readAlignment);
    alignmentClassifier_.classify(mateAlignment);
}

int SmallVariantAnalyzer::countReadsSupportingNode(graphtools::NodeId nodeId) const
{
    if (nodeId == ClassifierOfAlignmentsToVariant::kInvalidNodeId)
    {
        return alignmentClassifier_.numBypassingReads();
    }

    const auto& spanningCounts = alignmentClassifier_.countsOfSpanningReads();
    const auto& upstreamFlankingCounts = alignmentClassifier_.countsOfReadsFlankingUpstream();
    const auto& downstreamFlankingCounts = alignmentClassifier_.countsOfReadsFlankingDownstream();

    const int numReadsSupportingUpstreamFlank = upstreamFlankingCounts.countOf(nodeId) + spanningCounts.countOf(nodeId);
    const int numReadsSupportingDownstreamFlank
        = downstreamFlankingCounts.countOf(nodeId) + spanningCounts.countOf(nodeId);

    return (numReadsSupportingUpstreamFlank + numReadsSupportingDownstreamFlank) / 2;
}

std::unique_ptr<VariantFindings> SmallVariantAnalyzer::analyze() const
{
    NodeId refNode = optionalRefNode_ ? *optionalRefNode_ : ClassifierOfAlignmentsToVariant::kInvalidNodeId;
    NodeId altNode = ClassifierOfAlignmentsToVariant::kInvalidNodeId;

    switch (variantSubtype_)
    {
    case VariantSubtype::kInsertion:
        altNode = nodeIds_.front();
        break;
    case VariantSubtype::kDeletion:
        altNode = ClassifierOfAlignmentsToVariant::kInvalidNodeId;
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

    auto genotype = smallVariantGenotyper_.genotype(refNodeSupport, altNodeSupport);

    AllelePresenceStatus refAlleleStatus = allelePresenceChecker_.check(refNodeSupport, altNodeSupport);
    AllelePresenceStatus altAlleleStatus = allelePresenceChecker_.check(altNodeSupport, refNodeSupport);

    return std::unique_ptr<VariantFindings>(
        new SmallVariantFindings(refNodeSupport, altNodeSupport, refAlleleStatus, altAlleleStatus, genotype));
}

}
