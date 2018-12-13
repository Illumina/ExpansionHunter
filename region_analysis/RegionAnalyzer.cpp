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

#include "region_analysis/RegionAnalyzer.hh"

#include <cassert>
#include <iostream>

#include "thirdparty/spdlog/spdlog.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "alignment/AlignmentFilters.hh"
#include "alignment/AlignmentTweakers.hh"
#include "alignment/GraphAlignmentOperations.hh"
#include "region_analysis/RepeatAnalyzer.hh"
#include "region_analysis/SmallVariantAnalyzer.hh"

namespace ehunter
{

using boost::optional;
using graphtools::Operation;
using graphtools::OperationType;
using graphtools::splitStringByDelimiter;
using reads::LinearAlignmentStats;
using reads::Read;
using std::list;
using std::string;
using std::vector;

static const string encodeReadPair(const Read& read, const Read& mate)
{
    return read.read_id + ": " + read.sequence + "\n" + mate.read_id + ": " + read.sequence;
}

static string indentMultilineString(const string str, int32_t indentation_len)
{
    string indented_str;
    const vector<string> lines = splitStringByDelimiter(str, '\n');
    for (auto& line : lines)
    {
        if (!indented_str.empty())
        {
            indented_str += '\n';
        }
        indented_str += string(indentation_len, ' ') + line;
    }

    return indented_str;
}

static void outputAlignedRead(const Read& read, GraphAlignment alignment, std::ostream& out)
{
    const int32_t indentationSize = 2;
    const string spacer(indentationSize, ' ');
    out << spacer << "- name: " << read.read_id << std::endl;
    out << spacer << "  path: " << alignment.path() << std::endl;
    out << spacer << "  graph_cigar: " << alignment.generateCigar() << std::endl;
    out << spacer << "  alignment: |" << std::endl;
    const string alignmentEncoding = prettyPrint(alignment, read.sequence);
    out << indentMultilineString(alignmentEncoding, 3 * indentationSize) << std::endl;
}

RegionAnalyzer::RegionAnalyzer(
    const LocusSpecification& regionSpec, SampleParameters sampleParams, HeuristicParameters heuristicParams,
    std::ostream& alignmentStream)
    : regionSpec_(regionSpec)
    , sampleParams_(sampleParams)
    , heuristicParams_(heuristicParams)
    , alignmentStream_(alignmentStream)
    , orientationPredictor_(sampleParams.readLength(), &regionSpec_.regionGraph())
    , graphAligner_(
          &regionSpec_.regionGraph(), heuristicParams.alignerType(), heuristicParams_.kmerLenForAlignment(),
          heuristicParams_.paddingLength(), heuristicParams_.seedAffixTrimLength())
{
    verboseLogger_ = spdlog::get("verbose");

    for (const auto& variantSpec : regionSpec_.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kRepeat)
        {
            const auto& graph = regionSpec_.regionGraph();
            const int repeatNodeId = variantSpec.nodes().front();
            const string& repeatUnit = graph.nodeSeq(repeatNodeId);
            const int repeatUnitLength = repeatUnit.length();
            const int maxNumUnitsInRead = std::ceil(sampleParams_.readLength() / static_cast<double>(repeatUnitLength));

            weightedPurityCalculators.emplace(std::make_pair(repeatUnit, WeightedPurityCalculator(repeatUnit)));

            if (variantSpec.classification().subtype == VariantSubtype::kRareRepeat)
            {
                if (optionalUnitOfRareRepeat_)
                {
                    const string errorMessage
                        = "Region " + regionSpec_.regionId() + " is not permitted to have more than one rare variant";
                    throw std::logic_error(errorMessage);
                }
                optionalUnitOfRareRepeat_ = repeatUnit;
            }

            variantAnalyzerPtrs_.emplace_back(new RepeatAnalyzer(
                variantSpec.id(), regionSpec.expectedAlleleCount(), regionSpec.regionGraph(), repeatNodeId,
                sampleParams_.haplotypeDepth(), maxNumUnitsInRead));
        }
        else if (variantSpec.classification().type == VariantType::kSmallVariant)
        {
            variantAnalyzerPtrs_.emplace_back(new SmallVariantAnalyzer(
                variantSpec.id(), variantSpec.classification().subtype, regionSpec.expectedAlleleCount(),
                regionSpec.regionGraph(), variantSpec.nodes(), variantSpec.optionalRefNode(),
                sampleParams_.haplotypeDepth()));
        }
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw std::logic_error("Missing logic to create an analyzer for " + encoding.str());
        }
    }
}

void RegionAnalyzer::processMates(reads::Read read, reads::Read mate)
{
    optional<GraphAlignment> readAlignment = alignRead(read);
    optional<GraphAlignment> mateAlignment = alignRead(mate);

    int kMinNonRepeatAlignmentScore = sampleParams_.readLength() / 7.5;
    kMinNonRepeatAlignmentScore = std::max(kMinNonRepeatAlignmentScore, 3);
    if (!checkIfLocallyPlacedReadPair(readAlignment, mateAlignment, kMinNonRepeatAlignmentScore))
    {
        if (verboseLogger_)
        {
            verboseLogger_->info("Not locally placed read pair");
            verboseLogger_->info(encodeReadPair(read, mate));
        }
        return;
    }

    if (readAlignment && mateAlignment)
    {
        outputAlignedRead(read, *readAlignment, alignmentStream_);
        outputAlignedRead(mate, *mateAlignment, alignmentStream_);

        for (auto& variantAnalyzerPtr : variantAnalyzerPtrs_)
        {
            variantAnalyzerPtr->processMates(read, *readAlignment, mate, *mateAlignment);
        }
    }
    else if (verboseLogger_)
    {
        const string readStatus = readAlignment ? "Able" : "Unable";
        const string mateStatus = mateAlignment ? "Able" : "Unable";
        verboseLogger_->info("{} to align {}: {}", readStatus, read.readId(), read.sequence);
        verboseLogger_->info("{} to align {}: {}", mateStatus, mate.readId(), mate.sequence);
    }
}

void RegionAnalyzer::processOfftargetMates(reads::Read read1, reads::Read read2)
{
    if (!optionalUnitOfRareRepeat_)
    {
        const string errorMessage
            = "Cannot process offtarget mates for " + regionSpec_.regionId() + " because repeat unit is not set";
        throw std::logic_error(errorMessage);
    }

    const string& repeatUnit = *optionalUnitOfRareRepeat_;

    const auto& weightedPurityCalculator = weightedPurityCalculators.at(repeatUnit);
    const bool isFirstReadInrepeat = weightedPurityCalculator.score(read1.sequence) >= 0.90;
    const bool isSecondReadInrepeat = weightedPurityCalculator.score(read2.sequence) >= 0.90;

    if (isFirstReadInrepeat && isSecondReadInrepeat)
    {
        std::cerr << "Found IRR pair " << read1.fragmentId() << std::endl;
        processMates(std::move(read1), std::move(read2));
    }
}

boost::optional<GraphAlignment> RegionAnalyzer::alignRead(Read& read) const
{
    OrientationPrediction predictedOrientation = orientationPredictor_.predict(read.sequence);

    if (predictedOrientation == OrientationPrediction::kAlignsInReverseComplementOrientation)
    {
        read.sequence = graphtools::reverseComplement(read.sequence);
    }
    else if (predictedOrientation == OrientationPrediction::kDoesNotAlign)
    {
        return boost::optional<GraphAlignment>();
    }

    const list<GraphAlignment> alignments = graphAligner_.align(read.sequence);

    if (alignments.empty())
    {
        return boost::optional<GraphAlignment>();
    }

    GraphAlignment canonicalAlignment = computeCanonicalAlignment(alignments);

    if (checkIfPassesAlignmentFilters(canonicalAlignment))
    {
        // const int kShrinkLength = 10;
        // shrinkUncertainPrefix(kShrinkLength, read.sequence, canonicalAlignment);
        // shrinkUncertainSuffix(kShrinkLength, read.sequence, canonicalAlignment);

        return canonicalAlignment;
    }
    else
    {
        return boost::optional<GraphAlignment>();
    }
}

bool RegionAnalyzer::checkIfPassesAlignmentFilters(const GraphAlignment& alignment) const
{
    const Operation& firstOperation = alignment.alignments().front().operations().front();
    const int frontSoftclipLen = firstOperation.type() == OperationType::kSoftclip ? firstOperation.queryLength() : 0;

    const Operation& lastOperation = alignment.alignments().back().operations().back();
    const int backSoftclipLen = lastOperation.type() == OperationType::kSoftclip ? lastOperation.queryLength() : 0;

    const int clippedQueryLength = alignment.queryLength() - frontSoftclipLen - backSoftclipLen;
    const int referenceLength = alignment.referenceLength();

    const int percentQueryMatches = (100 * alignment.numMatches()) / clippedQueryLength;
    const int percentReferenceMatches = (100 * alignment.numMatches()) / referenceLength;

    if (percentQueryMatches >= 80 && percentReferenceMatches >= 80)
    {
        return true;
    }
    else
    {
        return false;
    }
}

RegionFindings RegionAnalyzer::genotype()
{
    RegionFindings regionResults;

    for (auto& variantAnalyzerPtr : variantAnalyzerPtrs_)
    {
        std::unique_ptr<VariantFindings> variantFindingsPtr = variantAnalyzerPtr->analyze();
        regionResults.emplace(std::make_pair(variantAnalyzerPtr->variantId(), std::move(variantFindingsPtr)));
    }

    return regionResults;
}

vector<std::unique_ptr<RegionAnalyzer>> initializeRegionAnalyzers(
    const RegionCatalog& RegionCatalog, const SampleParameters& sampleParams,
    const HeuristicParameters& heuristicParams, std::ostream& alignmentStream)
{
    vector<std::unique_ptr<RegionAnalyzer>> regionAnalyzers;

    for (const auto& regionIdAndRegionSpec : RegionCatalog)
    {
        const LocusSpecification& regionSpec = regionIdAndRegionSpec.second;

        // alignmentStream << regionSpec.regionId() << ":" << std::endl;
        regionAnalyzers.emplace_back(new RegionAnalyzer(regionSpec, sampleParams, heuristicParams, alignmentStream));
    }

    return regionAnalyzers;
}

}
