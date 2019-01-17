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

#include <sstream>

using boost::optional;
using graphtools::AlignmentWriter;
using graphtools::Operation;
using graphtools::OperationType;
using graphtools::splitStringByDelimiter;
using std::list;
using std::string;
using std::vector;

static const string encodeReadPair(const Read& read, const Read& mate)
{
    std::stringstream encoding;
    encoding << read.readId() << ": " << read.sequence() << "\n" << mate.readId() << ": " << read.sequence();
    return encoding.str();
}

RegionAnalyzer::RegionAnalyzer(
    const LocusSpecification& regionSpec, HeuristicParameters heuristicParams,
    graphtools::AlignmentWriter& alignmentWriter)
    : regionSpec_(regionSpec)
    , heuristicParams_(heuristicParams)
    , alignmentWriter_(alignmentWriter)
    , orientationPredictor_(&regionSpec_.regionGraph())
    , graphAligner_(
          &regionSpec_.regionGraph(), heuristicParams.alignerType(), heuristicParams_.kmerLenForAlignment(),
          heuristicParams_.paddingLength(), heuristicParams_.seedAffixTrimLength())
    , console_(spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console"))

{
    for (const auto& variantSpec : regionSpec_.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kRepeat)
        {
            const auto& graph = regionSpec_.regionGraph();
            const int repeatNodeId = variantSpec.nodes().front();
            const string& repeatUnit = graph.nodeSeq(repeatNodeId);

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
                variantSpec.id(), regionSpec.expectedAlleleCount(), regionSpec.regionGraph(), repeatNodeId));
        }
        else if (variantSpec.classification().type == VariantType::kSmallVariant)
        {
            variantAnalyzerPtrs_.emplace_back(new SmallVariantAnalyzer(
                variantSpec.id(), variantSpec.classification().subtype, regionSpec.expectedAlleleCount(),
                regionSpec.regionGraph(), variantSpec.nodes(), variantSpec.optionalRefNode()));
        }
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw std::logic_error("Missing logic to create an analyzer for " + encoding.str());
        }
    }
}

void RegionAnalyzer::processMates(Read read, Read mate)
{
    optional<GraphAlignment> readAlignment = alignRead(read);
    optional<GraphAlignment> mateAlignment = alignRead(mate);

    int kMinNonRepeatAlignmentScore = read.sequence().length() / 7.5;
    kMinNonRepeatAlignmentScore = std::max(kMinNonRepeatAlignmentScore, 10);
    if (!checkIfLocallyPlacedReadPair(readAlignment, mateAlignment, kMinNonRepeatAlignmentScore))
    {
        console_->debug("Not locally placed read pair\n" + encodeReadPair(read, mate));
        return;
    }

    if (readAlignment && mateAlignment)
    {
        alignmentWriter_.write(
            regionSpec_.regionId(), read.fragmentId(), read.sequence(), read.isFirstMate(), *readAlignment);
        alignmentWriter_.write(
            regionSpec_.regionId(), mate.fragmentId(), mate.sequence(), mate.isFirstMate(), *mateAlignment);

        for (auto& variantAnalyzerPtr : variantAnalyzerPtrs_)
        {
            variantAnalyzerPtr->processMates(read, *readAlignment, mate, *mateAlignment);
        }
    }
    else
    {
        const string readStatus = readAlignment ? "Able" : "Unable";
        const string mateStatus = mateAlignment ? "Able" : "Unable";
        console_->debug("{} to align {}: {}", readStatus, read.readId(), read.sequence());
        console_->debug("{} to align {}: {}", mateStatus, mate.readId(), mate.sequence());
    }
}

void RegionAnalyzer::processOfftargetMates(Read read1, Read read2)
{
    if (!optionalUnitOfRareRepeat_)
    {
        const string errorMessage
            = "Cannot process offtarget mates for " + regionSpec_.regionId() + " because repeat unit is not set";
        throw std::logic_error(errorMessage);
    }

    const string& repeatUnit = *optionalUnitOfRareRepeat_;

    const auto& weightedPurityCalculator = weightedPurityCalculators.at(repeatUnit);
    const bool isFirstReadInrepeat = weightedPurityCalculator.score(read1.sequence()) >= 0.90;
    const bool isSecondReadInrepeat = weightedPurityCalculator.score(read2.sequence()) >= 0.90;

    if (isFirstReadInrepeat && isSecondReadInrepeat)
    {
        processMates(std::move(read1), std::move(read2));
    }
}

boost::optional<GraphAlignment> RegionAnalyzer::alignRead(Read& read) const
{
    OrientationPrediction predictedOrientation = orientationPredictor_.predict(read.sequence());

    if (predictedOrientation == OrientationPrediction::kAlignsInReverseComplementOrientation)
    {
        read.setSequence(graphtools::reverseComplement(read.sequence()));
    }
    else if (predictedOrientation == OrientationPrediction::kDoesNotAlign)
    {
        return boost::optional<GraphAlignment>();
    }

    const list<GraphAlignment> alignments = graphAligner_.align(read.sequence());

    if (alignments.empty())
    {
        return boost::optional<GraphAlignment>();
    }

    return computeCanonicalAlignment(alignments);
}

RegionFindings RegionAnalyzer::analyze(const SampleParameters& params)
{
    RegionFindings regionResults;

    for (auto& variantAnalyzerPtr : variantAnalyzerPtrs_)
    {
        std::unique_ptr<VariantFindings> variantFindingsPtr = variantAnalyzerPtr->analyze(params);
        regionResults.emplace(std::make_pair(variantAnalyzerPtr->variantId(), std::move(variantFindingsPtr)));
    }

    return regionResults;
}

vector<std::unique_ptr<RegionAnalyzer>> initializeRegionAnalyzers(
    const RegionCatalog& RegionCatalog, const HeuristicParameters& heuristicParams, AlignmentWriter& bamletWriter)
{
    vector<std::unique_ptr<RegionAnalyzer>> regionAnalyzers;

    for (const auto& regionIdAndRegionSpec : RegionCatalog)
    {
        const LocusSpecification& regionSpec = regionIdAndRegionSpec.second;

        regionAnalyzers.emplace_back(new RegionAnalyzer(regionSpec, heuristicParams, bamletWriter));
    }

    return regionAnalyzers;
}

}
