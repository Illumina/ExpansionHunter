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

#include "region_analysis/RepeatAnalyzer.hh"

#include <iterator>

#include "thirdparty/spdlog/spdlog.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "alignment/OperationsOnAlignments.hh"
#include "genotyping/RepeatGenotyper.hh"

namespace ehunter
{

using boost::optional;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::prettyPrint;
using graphtools::splitStringByDelimiter;
using std::list;
using std::set;
using std::string;
using std::unique_ptr;
using std::vector;

void RepeatAnalyzer::processMates(
    const Read& read, const list<GraphAlignment>& readAlignments, const Read& mate,
    const list<GraphAlignment>& mateAlignments)
{
    ReadSummaryForStr strRead = alignmentClassifier_.classifyRead(read.sequence(), readAlignments);
    if (strRead.hasAlignments())
    {
        const StrAlignment& strAlignment = strRead.alignments().front();
        console_->trace("{} is {} for variant {}", read.readId(), strAlignment.type(), variantId_);
        summarizeAlignmentsToReadCounts(strAlignment);
    }
    else
    {
        console_->debug(
            "Could not confidently align {} to repeat node {} of {}", read.readId(), repeatNodeId(), variantId_);
    }

    ReadSummaryForStr strMate = alignmentClassifier_.classifyRead(mate.sequence(), mateAlignments);
    if (strMate.hasAlignments())
    {
        const StrAlignment& strAlignment = strMate.alignments().front();
        console_->trace("{} is {} for variant {}", mate.readId(), strAlignment.type(), variantId_);
        summarizeAlignmentsToReadCounts(strAlignment);
    }
    else
    {
        console_->debug(
            "Could not confidently align {} to repeat node {} of {}", mate.readId(), repeatNodeId(), variantId_);
    }
}

void RepeatAnalyzer::summarizeAlignmentsToReadCounts(const StrAlignment& strAlignment)
{
    switch (strAlignment.type())
    {
    case StrAlignment::Type::kSpanning:
        countsOfSpanningReads_.incrementCountOf(strAlignment.numUnits());
        break;
    case StrAlignment::Type::kFlanking:
        countsOfFlankingReads_.incrementCountOf(strAlignment.numUnits());
        break;
    case StrAlignment::Type::kInrepeat:
        countsOfInrepeatReads_.incrementCountOf(strAlignment.numUnits());
        break;
    default:
        break;
    }
}

static vector<int32_t> generateCandidateAlleleSizes(
    const CountTable& spanningTable, const CountTable& flankingTable, const CountTable& inrepeatTable)
{
    vector<int32_t> candidateSizes = spanningTable.getElementsWithNonzeroCounts();
    const int32_t longestSpanning
        = candidateSizes.empty() ? 0 : *std::max_element(candidateSizes.begin(), candidateSizes.end());

    const vector<int32_t> flankingSizes = flankingTable.getElementsWithNonzeroCounts();
    const int32_t longestFlanking
        = flankingSizes.empty() ? 0 : *std::max_element(flankingSizes.begin(), flankingSizes.end());

    const vector<int32_t> inrepeatSizes = inrepeatTable.getElementsWithNonzeroCounts();
    const int32_t longestInrepeat
        = inrepeatSizes.empty() ? 0 : *std::max_element(inrepeatSizes.begin(), inrepeatSizes.end());

    const int32_t longestNonSpanning = std::max(longestFlanking, longestInrepeat);

    if (longestSpanning < longestNonSpanning)
    {
        candidateSizes.push_back(longestNonSpanning);
    }

    return candidateSizes;
}

unique_ptr<VariantFindings> RepeatAnalyzer::analyze(const LocusStats& stats) const
{

    const int32_t repeatUnitLength = repeatUnit_.length();
    const double propCorrectMolecules = 0.97;
    const int maxNumUnitsInRead = std::ceil(stats.meanReadLength() / static_cast<double>(repeatUnitLength));

    auto truncatedSpanningTable = collapseTopElements(countsOfSpanningReads_, maxNumUnitsInRead);
    auto truncatedFlankingTable = collapseTopElements(countsOfFlankingReads_, maxNumUnitsInRead);
    auto truncatedInrepeatTable = collapseTopElements(countsOfInrepeatReads_, maxNumUnitsInRead);

    const vector<int32_t> candidateAlleleSizes
        = generateCandidateAlleleSizes(truncatedSpanningTable, truncatedFlankingTable, truncatedInrepeatTable);

    const double haplotypeDepth = stats.alleleCount() == AlleleCount::kTwo ? stats.depth() / 2 : stats.depth();

    RepeatGenotyper repeatGenotyper(
        haplotypeDepth, stats.alleleCount(), repeatUnitLength, maxNumUnitsInRead, propCorrectMolecules,
        truncatedSpanningTable, truncatedFlankingTable, truncatedInrepeatTable, countOfInrepeatReadPairs_);

    optional<RepeatGenotype> repeatGenotype = repeatGenotyper.genotypeRepeat(candidateAlleleSizes);

    unique_ptr<VariantFindings> variantFindingsPtr(
        new RepeatFindings(truncatedSpanningTable, truncatedFlankingTable, truncatedInrepeatTable, repeatGenotype));
    return variantFindingsPtr;
}
}
