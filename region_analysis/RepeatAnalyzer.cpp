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

#include "region_analysis/RepeatAnalyzer.hh"

#include "thirdparty/spdlog/spdlog.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "alignment/AlignmentFilters.hh"
#include "alignment/GraphAlignmentOperations.hh"
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
using std::vector;

static bool checkIfAlignmentIsConfident(
    NodeId repeatNodeId, const GraphAlignment& alignment, const RepeatAlignmentStats& alignmentStats)
{
    if (!checkIfPassesAlignmentFilters(alignment))
    {
        return false;
    }

    const bool doesReadAlignWellOverLeftFlank = checkIfUpstreamAlignmentIsGood(repeatNodeId, alignment);
    const bool doesReadAlignWellOverRightFlank = checkIfDownstreamAlignmentIsGood(repeatNodeId, alignment);

    if (alignmentStats.canonicalAlignmentType() == AlignmentType::kFlanksRepeat)
    {
        if (!doesReadAlignWellOverLeftFlank && !doesReadAlignWellOverRightFlank)
        {
            return false;
        }
    }

    if (alignmentStats.canonicalAlignmentType() == AlignmentType::kSpansRepeat)
    {
        if (!doesReadAlignWellOverRightFlank || !doesReadAlignWellOverRightFlank)
        {
            return false;
        }
    }

    return true;
}

void RepeatAnalyzer::processMates(
    const Read& read, const GraphAlignment& readAlignment, const Read& mate, const GraphAlignment& mateAlignment)
{
    RepeatAlignmentStats readAlignmentStats = classifyReadAlignment(readAlignment);
    RepeatAlignmentStats mateAlignmentStats = classifyReadAlignment(mateAlignment);

    const bool isReadAlignmentConfident
        = checkIfAlignmentIsConfident(repeatNodeId(), readAlignment, readAlignmentStats);
    const bool isMateAlignmentConfident
        = checkIfAlignmentIsConfident(repeatNodeId(), mateAlignment, mateAlignmentStats);

    if (isReadAlignmentConfident)
    {
        // std::cerr << read.readId() << " is " << readAlignmentStats.canonicalAlignmentType() << " for variant "
        //          << variantId_ << std::endl;
        summarizeAlignmentsToReadCounts(readAlignmentStats);
    }
    else
    {
        // std::cerr << read.readId() << " could not be confidently aligned " << std::endl;
        // std::cerr << prettyPrint(readAlignment, read.sequence()) << std::endl;
        console_->debug(
            "Not a confident alignment for repeat node {}\n{}", repeatNodeId(),
            prettyPrint(readAlignment, read.sequence()));
    }

    if (isMateAlignmentConfident)
    {
        // std::cerr << mate.readId() << " is " << mateAlignmentStats.canonicalAlignmentType() << " for variant "
        //          << std::endl;
        summarizeAlignmentsToReadCounts(mateAlignmentStats);
    }
    else
    {
        // std::cerr << mate.readId() << " could not be confidently aligned " << std::endl;
        // std::cerr << prettyPrint(mateAlignment, mate.sequence()) << std::endl;
        console_->debug(
            "Not a confident alignment for repeat node {}\n{}", repeatNodeId(),
            prettyPrint(mateAlignment, mate.sequence()));
    }
}

RepeatAlignmentStats RepeatAnalyzer::classifyReadAlignment(const graphtools::GraphAlignment& alignment)
{
    AlignmentType alignmentType = alignmentClassifier_.Classify(alignment);
    int numRepeatUnitsOverlapped = countFullOverlaps(repeatNodeId(), alignment);

    return RepeatAlignmentStats(alignment, alignmentType, numRepeatUnitsOverlapped);
}

void RepeatAnalyzer::summarizeAlignmentsToReadCounts(const RepeatAlignmentStats& repeatAlignmentStats)
{
    switch (repeatAlignmentStats.canonicalAlignmentType())
    {
    case AlignmentType::kSpansRepeat:
        countsOfSpanningReads_.incrementCountOf(repeatAlignmentStats.numRepeatUnitsSpanned());
        break;
    case AlignmentType::kFlanksRepeat:
        countsOfFlankingReads_.incrementCountOf(repeatAlignmentStats.numRepeatUnitsSpanned());
        break;
    case AlignmentType::kInsideRepeat:
        countsOfInrepeatReads_.incrementCountOf(repeatAlignmentStats.numRepeatUnitsSpanned());
        break;
    default:
        break;
    }
}

std::unique_ptr<VariantFindings> RepeatAnalyzer::analyze(const SampleParameters& params) const
{
    //     numRepeatUnitsOverlapped = std::min(numRepeatUnitsOverlapped, maxNumUnitsInRead_);

    const vector<int32_t> candidateAlleleSizes = generateCandidateAlleleSizes();

    const int32_t repeatUnitLength = repeatUnit_.length();
    const double propCorrectMolecules = 0.97;
    const int maxNumUnitsInRead = std::ceil(params.readLength() / static_cast<double>(repeatUnitLength));

    RepeatGenotyper repeatGenotyper(
        params.haplotypeDepth(), expectedAlleleCount_, repeatUnitLength, maxNumUnitsInRead, propCorrectMolecules,
        countsOfSpanningReads_, countsOfFlankingReads_, countsOfInrepeatReads_);

    optional<RepeatGenotype> repeatGenotype = repeatGenotyper.genotypeRepeat(candidateAlleleSizes);

    std::unique_ptr<VariantFindings> variantFiningsPtr(
        new RepeatFindings(countsOfSpanningReads_, countsOfFlankingReads_, countsOfInrepeatReads_, repeatGenotype));
    return variantFiningsPtr;
}

vector<int32_t> RepeatAnalyzer::generateCandidateAlleleSizes() const
{
    vector<int32_t> candidateAlleleSizes = countsOfSpanningReads_.getElementsWithNonzeroCounts();
    const int32_t longestSpanning = candidateAlleleSizes.empty()
        ? 0
        : *std::max_element(candidateAlleleSizes.begin(), candidateAlleleSizes.end());

    const vector<int32_t> repeatSizesInFlankingReads = countsOfFlankingReads_.getElementsWithNonzeroCounts();
    const int32_t longestFlanking = repeatSizesInFlankingReads.empty()
        ? 0
        : *std::max_element(repeatSizesInFlankingReads.begin(), repeatSizesInFlankingReads.end());

    const vector<int32_t> repeatSizesInInrepeatReads = countsOfInrepeatReads_.getElementsWithNonzeroCounts();
    const int32_t longestInrepeat = repeatSizesInInrepeatReads.empty()
        ? 0
        : *std::max_element(repeatSizesInInrepeatReads.begin(), repeatSizesInInrepeatReads.end());

    const int32_t longestNonSpanning = std::max(longestFlanking, longestInrepeat);

    if (longestSpanning < longestNonSpanning)
    {
        candidateAlleleSizes.push_back(longestNonSpanning);
    }

    return candidateAlleleSizes;
}

}
