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

#include "locus/SmallVariantAnalyzer.hh"

#include <vector>

#include <boost/optional.hpp>
#include <boost/smart_ptr/make_unique.hpp>

using boost::make_unique;
using boost::optional;
using graphtools::NodeId;
using std::string;
using std::vector;

namespace ehunter
{

void SmallVariantAnalyzer::processMates(
    const Read& /*read*/, const graphtools::GraphAlignment& readAlignment, const Read& /*mate*/,
    const graphtools::GraphAlignment& mateAlignment)
{
    alignmentStatsCalculator_.inspect(readAlignment);
    alignmentStatsCalculator_.inspect(mateAlignment);

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

std::unique_ptr<VariantFindings> SmallVariantAnalyzer::analyze(const LocusStats& stats)
{
    if (isLowDepth(stats))
    {
        auto refStatus = AlleleCheckSummary(AlleleStatus::kUncertain, 0);
        auto altStatus = AlleleCheckSummary(AlleleStatus::kUncertain, 0);
        return make_unique<SmallVariantFindings>(
            0, 0, refStatus, altStatus, stats.alleleCount(), boost::none, GenotypeFilter::kLowDepth);
    }

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

    int refNodeSupport = countReadsSupportingNode(refNode);
    int altNodeSupport = countReadsSupportingNode(altNode);

    const double haplotypeDepth = stats.alleleCount() == AlleleCount::kTwo ? stats.depth() / 2 : stats.depth();
    const int minBreakpointSpanningReads = stats.alleleCount() == AlleleCount::kTwo
        ? genotyperParams_.minBreakpointSpanningReads
        : (genotyperParams_.minBreakpointSpanningReads / 2);

    SmallVariantGenotyper smallVariantGenotyper(haplotypeDepth, stats.alleleCount());
    auto genotype = smallVariantGenotyper.genotype(refNodeSupport, altNodeSupport);

    auto refAlleleStatus = allelePresenceChecker_.check(haplotypeDepth, refNodeSupport, altNodeSupport);
    auto altAlleleStatus = allelePresenceChecker_.check(haplotypeDepth, altNodeSupport, refNodeSupport);

    GraphVariantAlignmentStats alignmentStats = alignmentStatsCalculator_.getStats();
    const bool insufficientBreakpointCoverage
        = alignmentStats.numReadsSpanningLeftBreakpoint() < minBreakpointSpanningReads
        || alignmentStats.numReadsSpanningRightBreakpoint() < minBreakpointSpanningReads;

    auto genotypeFilter = GenotypeFilter();
    if ((variantSubtype_ != VariantSubtype::kSMN) && insufficientBreakpointCoverage)
    {
        genotypeFilter = genotypeFilter | GenotypeFilter::kLowDepth;
    }

    return make_unique<SmallVariantFindings>(
        refNodeSupport, altNodeSupport, refAlleleStatus, altAlleleStatus, stats.alleleCount(), genotype,
        genotypeFilter);
}

}
