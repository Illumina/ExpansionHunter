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

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "alignment/AlignmentFilters.hh"
#include "alignment/GraphAlignmentOperations.hh"
#include "genotyping/RepeatGenotyper.hh"

using boost::optional;
using graphtools::GraphAlignment;
using graphtools::splitStringByDelimiter;
using reads::RepeatAlignmentStats;
using std::list;
using std::string;
using std::vector;

#define INTENTIONALLY_UNUSED(x) (void)(x)

void RepeatAnalyzer::processMates(
    const string& readName, const std::string& readSequence, boost::optional<GraphAlignment> readAlignment,
    const std::string& mateName, const std::string& mateSequence, boost::optional<GraphAlignment> mateAlignment)
{
    INTENTIONALLY_UNUSED(readName);
    INTENTIONALLY_UNUSED(readSequence);
    INTENTIONALLY_UNUSED(mateName);
    INTENTIONALLY_UNUSED(mateSequence);

    const int readLength = maxNumUnitsInRead_ * repeatUnit_.length();
    int kMinNonRepeatAlignmentScore = readLength / 7.5;
    kMinNonRepeatAlignmentScore = std::max(kMinNonRepeatAlignmentScore, 3);
    if (!checkIfLocallyPlacedReadPair(readAlignment, mateAlignment, kMinNonRepeatAlignmentScore))
    {
        return;
    }

    if (readAlignment)
    {
        RepeatAlignmentStats readAlignmentStats = classifyReadAlignment(*readAlignment);

        const bool leftFlankAlignmentIsGood = checkIfUpstreamAlignmentIsGood(repeatNodeId_, *readAlignment);
        const bool rightFlankAlignmentIsGood = checkIfDownstreamAlignmentIsGood(repeatNodeId_, *readAlignment);

        if (readAlignmentStats.canonicalAlignmentType() == AlignmentType::kFlanksRepeat)
        {
            if (!leftFlankAlignmentIsGood && !rightFlankAlignmentIsGood)
            {
                return;
            }
        }

        if (readAlignmentStats.canonicalAlignmentType() == AlignmentType::kSpansRepeat)
        {
            if (!leftFlankAlignmentIsGood || !rightFlankAlignmentIsGood)
            {
                return;
            }
        }

        summarizeAlignmentsToReadCounts(readAlignmentStats);
    }

    if (mateAlignment)
    {
        RepeatAlignmentStats mateAlignmentStats = classifyReadAlignment(*mateAlignment);

        const bool leftFlankAlignmentIsGood = checkIfUpstreamAlignmentIsGood(repeatNodeId_, *mateAlignment);
        const bool rightFlankAlignmentIsGood = checkIfDownstreamAlignmentIsGood(repeatNodeId_, *mateAlignment);

        if (mateAlignmentStats.canonicalAlignmentType() == AlignmentType::kFlanksRepeat)
        {
            if (!leftFlankAlignmentIsGood && !rightFlankAlignmentIsGood)
            {
                return;
            }
        }

        if (mateAlignmentStats.canonicalAlignmentType() == AlignmentType::kSpansRepeat)
        {
            if (!leftFlankAlignmentIsGood || !rightFlankAlignmentIsGood)
            {
                return;
            }
        }

        summarizeAlignmentsToReadCounts(mateAlignmentStats);
    }
}

RepeatAlignmentStats RepeatAnalyzer::classifyReadAlignment(const graphtools::GraphAlignment& alignment)
{
    AlignmentType alignmentType = alignmentClassifier_.Classify(alignment);
    int32_t numRepeatUnitsOverlapped = countFullOverlaps(repeatNodeId_, alignment);
    numRepeatUnitsOverlapped = std::min(numRepeatUnitsOverlapped, maxNumUnitsInRead_);

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
        countsOfRepeatReads_.incrementCountOf(repeatAlignmentStats.numRepeatUnitsSpanned());
        break;
    default:
        break;
    }
}

RepeatFindings RepeatAnalyzer::analyzeCommonRepeat() const
{
    const vector<int32_t> candidateAlleleSizes = generateCandidateAlleleSizes();
    optional<RepeatGenotype> repeatGenotype = genotypeCommonRepeat(candidateAlleleSizes);
    std::cerr << repeatId_ << " spanning = " << countsOfSpanningReads_ << " flanking = " << countsOfFlankingReads_
              << " inrepeat = " << countsOfRepeatReads_ << std::endl;
    return RepeatFindings(countsOfFlankingReads_, countsOfSpanningReads_, repeatGenotype);
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

    if (longestSpanning < longestFlanking)
    {
        candidateAlleleSizes.push_back(longestFlanking);
    }

    return candidateAlleleSizes;
}

optional<RepeatGenotype> RepeatAnalyzer::genotypeCommonRepeat(const vector<int32_t>& candidateAlleleSizes) const
{
    const int32_t repeatUnitLen = repeatUnit_.length();
    const double propCorrectMolecules = 0.97;

    RepeatGenotyper repeatGenotyper(
        haplotypeDepth_, expectedAlleleCount_, repeatUnitLen, maxNumUnitsInRead_, propCorrectMolecules,
        countsOfSpanningReads_, countsOfFlankingReads_, countsOfRepeatReads_);

    return repeatGenotyper.genotypeRepeat(candidateAlleleSizes);
}

bool RepeatAnalyzer::operator==(const RepeatAnalyzer& other) const
{
    return repeatId_ == other.repeatId_ && expectedAlleleCount_ == other.expectedAlleleCount_
        && repeatNodeId_ == other.repeatNodeId_ && maxNumUnitsInRead_ == other.maxNumUnitsInRead_
        && haplotypeDepth_ == other.haplotypeDepth_ && repeatUnit_ == other.repeatUnit_
        && alignmentClassifier_ == other.alignmentClassifier_ && countsOfFlankingReads_ == other.countsOfFlankingReads_
        && countsOfSpanningReads_ == other.countsOfSpanningReads_;
}
