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

#include "region/StrAnalyzer.hh"
#include "common/CountTable.hh"
#include "genotyping/RepeatGenotyper.hh"

using boost::optional;
using std::string;
using std::vector;

namespace ehunter
{

using ReadSummaries = vector<ReadSummaryForStr>;

static void populateCountTables(
    const ReadSummaries& readSummaries, CountTable& spanningReads, CountTable& flankingReads, CountTable& inrepeatReads)
{
    for (const auto& summary : readSummaries)
    {
        const auto& alignment = summary.alignments().front();

        switch (alignment.type())
        {
        case StrAlignment::Type::kSpanning:
            spanningReads.incrementCountOf(alignment.numUnits());
            break;
        case StrAlignment::Type::kFlanking:
            flankingReads.incrementCountOf(alignment.numUnits());
            break;
        case StrAlignment::Type::kInrepeat:
            inrepeatReads.incrementCountOf(alignment.numUnits());
            break;
        default:
            break;
        }
    }
}

static vector<int> generateCandidateAlleleSizes(
    const CountTable& spanningTable, const CountTable& flankingTable, const CountTable& inrepeatTable)
{
    vector<int> candidateSizes = spanningTable.getElementsWithNonzeroCounts();
    const int longestSpanning
        = candidateSizes.empty() ? 0 : *std::max_element(candidateSizes.begin(), candidateSizes.end());

    const vector<int> flankingSizes = flankingTable.getElementsWithNonzeroCounts();
    const int longestFlanking
        = flankingSizes.empty() ? 0 : *std::max_element(flankingSizes.begin(), flankingSizes.end());

    const vector<int> inrepeatSizes = inrepeatTable.getElementsWithNonzeroCounts();
    const int longestInrepeat
        = inrepeatSizes.empty() ? 0 : *std::max_element(inrepeatSizes.begin(), inrepeatSizes.end());

    const int longestNonSpanning = std::max(longestFlanking, longestInrepeat);

    if (longestSpanning < longestNonSpanning)
    {
        candidateSizes.push_back(longestNonSpanning);
    }

    return candidateSizes;
}

StrAnalyzer::StrAnalyzer(string variantId)
    : VariantAnalyzer(std::move(variantId))
{
}

std::unique_ptr<VariantFindings> StrAnalyzer::analyze(const LocusStats& stats) const
{
    CountTable spanningReads;
    CountTable flankingReads;
    CountTable inrepeatReads;
    populateCountTables(graphFeaturePtr_->readSummaries(), spanningReads, flankingReads, inrepeatReads);

    assert(graphFeaturePtr_);
    const string& motif = graphFeaturePtr_->motif();
    const int maxNumUnitsInRead = std::ceil(stats.meanReadLength() / static_cast<double>(motif.length()));
    auto truncatedSpanningTable = collapseTopElements(spanningReads, maxNumUnitsInRead);
    auto truncatedFlankingTable = collapseTopElements(flankingReads, maxNumUnitsInRead);
    auto truncatedInrepeatTable = collapseTopElements(inrepeatReads, maxNumUnitsInRead);

    const vector<int> candidateAlleleSizes
        = generateCandidateAlleleSizes(truncatedSpanningTable, truncatedFlankingTable, truncatedInrepeatTable);

    const double propCorrectMolecules = 0.97;
    const double haplotypeDepth = stats.alleleCount() == AlleleCount::kTwo ? stats.depth() / 2 : stats.depth();

    // TODO: Set correct count
    int countOfInrepeatReadPairs = 0;

    RepeatGenotyper repeatGenotyper(
        haplotypeDepth, stats.alleleCount(), motif.length(), maxNumUnitsInRead, propCorrectMolecules,
        truncatedSpanningTable, truncatedFlankingTable, truncatedInrepeatTable, countOfInrepeatReadPairs);

    optional<RepeatGenotype> genotype = repeatGenotyper.genotypeRepeat(candidateAlleleSizes);

    std::unique_ptr<VariantFindings> findingsPtr(
        new StrFindings(truncatedSpanningTable, truncatedFlankingTable, truncatedInrepeatTable, genotype));
    return findingsPtr;
}

}
