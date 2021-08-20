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
    optional<SmallVariantGenotype> genotype = boost::none;
    auto genotypeFilter = GenotypeFilter();

    int refNodeSupport = 0;
    int altNodeSupport = 0;

    AlleleCheckSummary refAlleleStatus(AlleleStatus::kUncertain, 0);
    AlleleCheckSummary altAlleleStatus(AlleleStatus::kUncertain, 0);

    if (stats.meanReadLength() == 0 || stats.depth() < genotyperParams_.minLocusCoverage)
    {
        genotypeFilter = genotypeFilter | GenotypeFilter::kLowDepth;
    }
    else
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

        refNodeSupport = countReadsSupportingNode(refNode);
        altNodeSupport = countReadsSupportingNode(altNode);

        const double haplotypeDepth = stats.alleleCount() == AlleleCount::kTwo ? stats.depth() / 2 : stats.depth();
        const int minBreakpointSpanningReads = stats.alleleCount() == AlleleCount::kTwo
            ? genotyperParams_.minBreakpointSpanningReads
            : (genotyperParams_.minBreakpointSpanningReads / 2);

        SmallVariantGenotyper smallVariantGenotyper(haplotypeDepth, stats.alleleCount());
        genotype = smallVariantGenotyper.genotype(refNodeSupport, altNodeSupport);

        refAlleleStatus = allelePresenceChecker_.check(haplotypeDepth, refNodeSupport, altNodeSupport);
        altAlleleStatus = allelePresenceChecker_.check(haplotypeDepth, altNodeSupport, refNodeSupport);

        GraphVariantAlignmentStats alignmentStats = alignmentStatsCalculator_.getStats();
        const bool insufficientBreakpointCoverage
            = alignmentStats.numReadsSpanningLeftBreakpoint() < minBreakpointSpanningReads
            || alignmentStats.numReadsSpanningRightBreakpoint() < minBreakpointSpanningReads;
        if ((variantSubtype_ != VariantSubtype::kSMN) && insufficientBreakpointCoverage)
        {
            genotypeFilter = genotypeFilter | GenotypeFilter::kLowDepth;
        }
    }

    return std::unique_ptr<VariantFindings>(new SmallVariantFindings(
        refNodeSupport, altNodeSupport, refAlleleStatus, altAlleleStatus, stats.alleleCount(), genotype,
        genotypeFilter));
}

}
