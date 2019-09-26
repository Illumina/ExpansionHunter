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

#include "workflow/GraphSmallVariantAnalyzer.hh"

#include "workflow/GraphSmallVariant.hh"

using graphtools::NodeId;
using std::shared_ptr;
using std::unique_ptr;
using std::vector;

namespace ehunter
{

GraphSmallVariantAnalyzer::GraphSmallVariantAnalyzer(
    std::shared_ptr<GraphSmallVariant> smallVariantFeature, std::string variantId, VariantSubtype variantSubtype,
    boost::optional<graphtools::NodeId> optionalRefNode)
    : GraphVariantAnalyzer(std::move(variantId))
    , smallVariantFeature_(std::move(smallVariantFeature))
    , variantSubtype_(variantSubtype)
    , optionalRefNode_(optionalRefNode)
    , allelePresenceChecker_(genotyperParams_.errorRate, genotyperParams_.likelihoodRatioThreshold)
{
}

unique_ptr<VariantFindings> GraphSmallVariantAnalyzer::analyze(const LocusStats& stats) const
{
    NodeId refNode = optionalRefNode_ ? *optionalRefNode_ : SmallVariantAlignmentClassifier::kInvalidNodeId;
    NodeId altNode = SmallVariantAlignmentClassifier::kInvalidNodeId;

    auto nodeIds = smallVariantFeature_->nodeIds();

    switch (variantSubtype_)
    {
    case VariantSubtype::kInsertion:
        altNode = nodeIds.front();
        break;
    case VariantSubtype::kDeletion:
        altNode = SmallVariantAlignmentClassifier::kInvalidNodeId;
        break;
    case VariantSubtype::kSwap:
        altNode = (refNode == nodeIds.front()) ? nodeIds.back() : nodeIds.front();
        break;
    case VariantSubtype::kSMN:
        if (refNode != nodeIds.front())
            throw std::logic_error("Invalid SMN specification");
        altNode = nodeIds.back();
        break;
    default:
        std::ostringstream encoding;
        encoding << variantSubtype_;
        throw std::logic_error("Invalid small variant subtype: " + encoding.str());
    }

    const int refNodeSupport = smallVariantFeature_->countReadsSupportingNode(refNode);
    const int altNodeSupport = smallVariantFeature_->countReadsSupportingNode(altNode);

    const double haplotypeDepth = stats.alleleCount() == AlleleCount::kTwo ? stats.depth() / 2 : stats.depth();

    SmallVariantGenotyper smallVariantGenotyper(haplotypeDepth, stats.alleleCount());
    auto genotype = smallVariantGenotyper.genotype(refNodeSupport, altNodeSupport);

    auto refAlleleStatus = allelePresenceChecker_.check(haplotypeDepth, refNodeSupport, altNodeSupport);
    auto altAlleleStatus = allelePresenceChecker_.check(haplotypeDepth, altNodeSupport, refNodeSupport);

    return std::unique_ptr<VariantFindings>(
        new SmallVariantFindings(refNodeSupport, altNodeSupport, refAlleleStatus, altAlleleStatus, genotype));
}

vector<shared_ptr<RegionModelFeature>> GraphSmallVariantAnalyzer::features() { return { smallVariantFeature_ }; }

}
