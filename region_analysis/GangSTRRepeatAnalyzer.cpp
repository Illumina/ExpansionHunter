//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Nima Mousavi <mousavi@ucsd.edu>
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

#include "region_analysis/GangSTRRepeatAnalyzer.hh"

#include "thirdparty/spdlog/spdlog.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "alignment/AlignmentFilters.hh"
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
    using std::string;
    using std::unique_ptr;
    using std::vector;

    static bool checkIfAlignmentIsConfident(
            NodeId repeatNodeId, const GraphAlignment& alignment, const GangSTRAlignmentStats& alignmentStats)
    {
        if (!checkIfPassesAlignmentFilters(alignment))
        {
            return false;
        }

        const bool doesReadAlignWellOverLeftFlank = checkIfUpstreamAlignmentIsGood(repeatNodeId, alignment);
        const bool doesReadAlignWellOverRightFlank = checkIfDownstreamAlignmentIsGood(repeatNodeId, alignment);

        if (alignmentStats.canonicalAlignmentType() == GangSTRAlignmentType::kFlanksLeft ||
                alignmentStats.canonicalAlignmentType() == GangSTRAlignmentType::kFlanksRight)
        {
            if (!doesReadAlignWellOverLeftFlank && !doesReadAlignWellOverRightFlank)
            {
                return false;
            }
        }

        if (alignmentStats.canonicalAlignmentType() == GangSTRAlignmentType::kSpansRepeat)
        {
            if (!doesReadAlignWellOverLeftFlank || !doesReadAlignWellOverRightFlank)
            {
                return false;
            }
        }

        return true;
    }

    void GangSTRRepeatAnalyzer::processMates(
            const Read& read, const GraphAlignment& readAlignment, const Read& mate, const GraphAlignment& mateAlignment)
    {
        // TODO Nima: combine these two function calls?
        GangSTRAlignmentStats readAlignmentStats = classifyReadAlignment(readAlignment, mateAlignment);
        GangSTRAlignmentStats mateAlignmentStats = classifyReadAlignment(mateAlignment, readAlignment);

        const bool isReadAlignmentConfident
                = checkIfAlignmentIsConfident(repeatNodeId(), readAlignment, readAlignmentStats);
        const bool isMateAlignmentConfident
                = checkIfAlignmentIsConfident(repeatNodeId(), mateAlignment, mateAlignmentStats);

        if (isReadAlignmentConfident)
        {
            console_->trace(
                    "{} is {} for variant {}", read.readId(), readAlignmentStats.canonicalAlignmentType(), variantId_);
            summarizeAlignmentsToReadCounts(readAlignmentStats);
        }
        else
        {
            console_->debug(
                    "Could not confidently aligned {} to repeat node {} of {}\n{}", read.readId(), repeatNodeId(), variantId_,
                    prettyPrint(readAlignment, read.sequence()));
        }

        if (isMateAlignmentConfident)
        {
            console_->trace(
                    "{} is {} for variant {}", mate.readId(), mateAlignmentStats.canonicalAlignmentType(), variantId_);
            summarizeAlignmentsToReadCounts(mateAlignmentStats);
        }
        else
        {
            console_->debug(
                    "Could not confidently align {} to repeat node {} of {}\n{}", mate.readId(), repeatNodeId(), variantId_,
                    prettyPrint(mateAlignment, mate.sequence()));
        }
    }

    GangSTRAlignmentStats GangSTRRepeatAnalyzer::classifyReadAlignment(const graphtools::GraphAlignment& alignment,
                                                                       const graphtools::GraphAlignment& mate)
    {
        GangSTRAlignmentType alignmentType = alignmentClassifier_.Classify(alignment);
        GangSTRAlignmentType mateType = alignmentClassifier_.Classify(mate);
        int numRepeatUnitsOverlapped = countFullOverlaps(repeatNodeId(), alignment);
        // TODO Nima: Calculate fragmentLength and mateDistanceToRepeat
        int32_t fragmentLength = -1;
        int32_t mateDistanceToRepeat = -1;
        return GangSTRAlignmentStats(alignment, alignmentType, mateType, numRepeatUnitsOverlapped, fragmentLength, mateDistanceToRepeat);
    }

    void GangSTRRepeatAnalyzer::summarizeAlignmentsToReadCounts(const GangSTRAlignmentStats& gangSTRAlignmentStats)
    {
        switch (gangSTRAlignmentStats.canonicalAlignmentType())
        {
            case GangSTRAlignmentType::kSpansRepeat:
                countsOfSpanningReads_.incrementCountOf(gangSTRAlignmentStats.numRepeatUnitsSpanned());
                break;
            case GangSTRAlignmentType::kFlanksLeft:
                countsOfFlankingReads_.incrementCountOf(gangSTRAlignmentStats.numRepeatUnitsSpanned());
                if (gangSTRAlignmentStats.mateAlignmentType() == GangSTRAlignmentType::kFlanksRight ||
                    gangSTRAlignmentStats.mateAlignmentType() == GangSTRAlignmentType::kRightOfRepeat) {
                    if (gangSTRAlignmentStats.fragmentLength() != -1) {
                        distanceOfTraversingPairs_.addElement(gangSTRAlignmentStats.fragmentLength());
                    }
                }
                break;
            case GangSTRAlignmentType::kFlanksRight:
                countsOfFlankingReads_.incrementCountOf(gangSTRAlignmentStats.numRepeatUnitsSpanned());
                if (gangSTRAlignmentStats.mateAlignmentType() == GangSTRAlignmentType::kFlanksLeft ||
                    gangSTRAlignmentStats.mateAlignmentType() == GangSTRAlignmentType::kLeftOfRepeat) {
                    if (gangSTRAlignmentStats.fragmentLength() != -1) {
                        distanceOfTraversingPairs_.addElement(gangSTRAlignmentStats.fragmentLength());
                    }
                }
                break;
            case GangSTRAlignmentType::kInsideRepeat:
                countsOfInrepeatReads_.incrementCountOf(gangSTRAlignmentStats.numRepeatUnitsSpanned());
                if (gangSTRAlignmentStats.mateAlignmentType() != GangSTRAlignmentType::kInsideRepeat) {
                    if (gangSTRAlignmentStats.mateDistanceToRepeat() != -1) {
                        distanceOfInrepeatMates_.addElement(gangSTRAlignmentStats.mateDistanceToRepeat());
                    }
                }
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

    unique_ptr<VariantFindings> GangSTRRepeatAnalyzer::analyze(const LocusStats& stats) const
    {

        const int32_t repeatUnitLength = repeatUnit_.length();
        const double propCorrectMolecules = 0.97;
        const int maxNumUnitsInRead = std::ceil(stats.meanReadLength() / static_cast<double>(repeatUnitLength));

        auto truncatedSpanningTable = collapseTopElements(countsOfSpanningReads_, maxNumUnitsInRead);
        auto truncatedFlankingTable = collapseTopElements(countsOfFlankingReads_, maxNumUnitsInRead);
        auto truncatedInrepeatTable = collapseTopElements(countsOfInrepeatReads_, maxNumUnitsInRead);

        const vector<int32_t> candidateAlleleSizes
                = generateCandidateAlleleSizes(truncatedSpanningTable, truncatedFlankingTable, truncatedInrepeatTable);

        const double haplotypeDepth = expectedAlleleCount_ == AlleleCount::kTwo ? stats.depth() / 2 : stats.depth();

        RepeatGenotyper repeatGenotyper(
                haplotypeDepth, expectedAlleleCount_, repeatUnitLength, maxNumUnitsInRead, propCorrectMolecules,
                truncatedSpanningTable, truncatedFlankingTable, truncatedInrepeatTable);

        optional<RepeatGenotype> repeatGenotype = repeatGenotyper.genotypeRepeat(candidateAlleleSizes);

        unique_ptr<VariantFindings> variantFindingsPtr(
                new RepeatFindings(truncatedSpanningTable, truncatedFlankingTable, truncatedInrepeatTable, repeatGenotype));
        return variantFindingsPtr;
    }

}
