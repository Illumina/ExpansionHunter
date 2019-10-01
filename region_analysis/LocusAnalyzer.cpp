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

#include "region_analysis/LocusAnalyzer.hh"

#include <cassert>
#include <iostream>

#include "thirdparty/spdlog/spdlog.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphutils/SequenceOperations.hh"

#include "alignment/AlignmentFilters.hh"
#include "alignment/AlignmentTweakers.hh"
#include "alignment/OperationsOnAlignments.hh"
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
using std::unique_ptr;
using std::vector;

/*
static const string encodeReadPair(const Read& read, const Read& mate)
{
    std::stringstream encoding;
    encoding << read.readId() << ": " << read.sequence() << "\n" << mate.readId() << ": " << read.sequence();
    return encoding.str();
}
*/

LocusAnalyzer::LocusAnalyzer(const LocusSpecification& locusSpec, graphtools::AlignmentWriter& alignmentWriter)
    : locusSpec_(locusSpec)
    , alignmentWriter_(alignmentWriter)
    , orientationPredictor_(&locusSpec_.regionGraph())
    , graphAligner_(
          &locusSpec_.regionGraph(), workflowContext_.heuristics().alignerType(),
          workflowContext_.heuristics().kmerLenForAlignment(), workflowContext_.heuristics().paddingLength(),
          workflowContext_.heuristics().seedAffixTrimLength())
    , statsCalculator_(locusSpec_.typeOfChromLocusLocatedOn(), locusSpec_.regionGraph())
    , console_(spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console"))

{
    for (const auto& variantSpec : locusSpec_.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kRepeat)
        {
            const auto& graph = locusSpec_.regionGraph();
            const int repeatNodeId = variantSpec.nodes().front();
            const string& repeatUnit = graph.nodeSeq(repeatNodeId);

            weightedPurityCalculators.emplace(std::make_pair(repeatUnit, WeightedPurityCalculator(repeatUnit)));

            if (variantSpec.classification().subtype == VariantSubtype::kRareRepeat)
            {
                if (optionalUnitOfRareRepeat_)
                {
                    const string errorMessage
                        = "Region " + locusSpec_.locusId() + " is not permitted to have more than one rare variant";
                    throw std::logic_error(errorMessage);
                }
                optionalUnitOfRareRepeat_ = repeatUnit;
            }

            variantAnalyzerPtrs_.emplace_back(
                new RepeatAnalyzer(variantSpec.id(), locusSpec.regionGraph(), repeatNodeId));
        }
        else if (variantSpec.classification().type == VariantType::kSmallVariant)
        {
            variantAnalyzerPtrs_.emplace_back(new SmallVariantAnalyzer(
                variantSpec.id(), variantSpec.classification().subtype, locusSpec.regionGraph(), variantSpec.nodes(),
                variantSpec.optionalRefNode(), locusSpec.genotyperParameters()));
        }
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw std::logic_error("Missing logic to create an analyzer for " + encoding.str());
        }
    }
}

void LocusAnalyzer::processMates(Read read, optional<Read> mate, RegionType regionType)
{
    if (regionType == RegionType::kTarget)
    {
        processOntargetMates(std::move(read), std::move(mate));
    }
    else if (mate)
    {
        processOfftargetMates(std::move(read), std::move(*mate));
    }
}

void LocusAnalyzer::processOntargetMates(Read read, optional<Read> mate)
{
    list<GraphAlignment> readAlignments = alignRead(read);
    list<GraphAlignment> mateAlignments;
    if (mate)
    {
        mateAlignments = alignRead(*mate);
    }

    int numMatchingBases = static_cast<int>(read.sequence().length() / 7.5);
    numMatchingBases = std::max(numMatchingBases, 10);
    LinearAlignmentParameters parameters;
    const int kMinNonRepeatAlignmentScore = numMatchingBases * parameters.matchScore;

    if (!checkIfComesFromGraphLocus(readAlignments, mateAlignments, kMinNonRepeatAlignmentScore))
    {
        if (mate && optionalUnitOfRareRepeat_)
        {
            processOfftargetMates(std::move(read), std::move(*mate));
        }

        return;
    }

    if (!readAlignments.empty())
    {
        statsCalculator_.inspect(readAlignments.front());
    }
    if (!mateAlignments.empty())
    {
        statsCalculator_.inspect(mateAlignments.front());
    }

    if (!readAlignments.empty() && !mateAlignments.empty())
    {
        runVariantAnalysis(read, readAlignments, *mate, mateAlignments);
    }
}

void LocusAnalyzer::processOfftargetMates(Read read, Read mate)
{
    if (!optionalUnitOfRareRepeat_)
    {
        const string errorMessage
            = "Cannot process offtarget mates for " + locusSpec_.locusId() + " because repeat unit is not set";
        throw std::logic_error(errorMessage);
    }

    const string& repeatUnit = *optionalUnitOfRareRepeat_;

    const auto& weightedPurityCalculator = weightedPurityCalculators.at(repeatUnit);
    const bool isFirstReadInrepeat = weightedPurityCalculator.score(read.sequence()) >= 0.90;
    const bool isSecondReadInrepeat = weightedPurityCalculator.score(mate.sequence()) >= 0.90;

    if (isFirstReadInrepeat && isSecondReadInrepeat)
    {
        int numAnalyzersFound = 0;
        for (auto& variantAnalyzerPtr : variantAnalyzerPtrs_)
        {
            auto repeatAnalyzerPtr = dynamic_cast<RepeatAnalyzer*>(variantAnalyzerPtr.get());
            if (repeatAnalyzerPtr != nullptr && repeatAnalyzerPtr->repeatUnit() == repeatUnit)
            {
                numAnalyzersFound++;
                repeatAnalyzerPtr->addInrepeatReadPair();
            }
        }

        if (numAnalyzersFound != 1)
        {
            const string errorMessage = "Encountered inconsistently-specified locus " + locusSpec_.locusId();
            throw std::logic_error(errorMessage);
        }
    }
}

void LocusAnalyzer::runVariantAnalysis(
    const Read& read, const GraphAlignments& readAlignments, const Read& mate, const GraphAlignments& mateAlignments)
{
    alignmentWriter_.write(
        locusSpec_.locusId(), read.fragmentId(), read.sequence(), read.isFirstMate(), read.isReversed(),
        mate.isReversed(), readAlignments.front());
    alignmentWriter_.write(
        locusSpec_.locusId(), mate.fragmentId(), mate.sequence(), mate.isFirstMate(), mate.isReversed(),
        read.isReversed(), mateAlignments.front());

    for (auto& variantAnalyzerPtr : variantAnalyzerPtrs_)
    {
        variantAnalyzerPtr->processMates(read, readAlignments, mate, mateAlignments);
    }
}

list<GraphAlignment> LocusAnalyzer::alignRead(Read& read) const
{
    OrientationPrediction predictedOrientation = orientationPredictor_.predict(read.sequence());

    if (predictedOrientation == OrientationPrediction::kAlignsInReverseComplementOrientation)
    {
        read.reverseComplement();
    }
    else if (predictedOrientation == OrientationPrediction::kDoesNotAlign)
    {
        return {};
    }

    return graphAligner_.align(read.sequence());
}

LocusFindings LocusAnalyzer::analyze(Sex sampleSex)
{
    LocusFindings locusFindings;

    locusFindings.optionalStats = statsCalculator_.estimate(sampleSex);
    if (locusFindings.optionalStats
        && locusFindings.optionalStats->depth() >= locusSpec().genotyperParameters().minLocusCoverage)
    {
        for (auto& variantAnalyzerPtr : variantAnalyzerPtrs_)
        {
            const LocusStats& locusStats = *locusFindings.optionalStats;
            std::unique_ptr<VariantFindings> variantFindingsPtr = variantAnalyzerPtr->analyze(locusStats);
            const string& variantId = variantAnalyzerPtr->variantId();
            locusFindings.findingsForEachVariant.emplace(variantId, std::move(variantFindingsPtr));
        }
    }
    else
    {
        console_->warn("Skipping locus {} due to low coverage", locusId());
    }

    return locusFindings;
}

vector<std::unique_ptr<LocusAnalyzer>>
initializeLocusAnalyzers(const RegionCatalog& RegionCatalog, AlignmentWriter& bamletWriter)
{
    vector<std::unique_ptr<LocusAnalyzer>> locusAnalyzers;

    for (const auto& locusIdAndRegionSpec : RegionCatalog)
    {
        const LocusSpecification& locusSpec = locusIdAndRegionSpec.second;

        locusAnalyzers.emplace_back(new LocusAnalyzer(locusSpec, bamletWriter));
    }

    return locusAnalyzers;
}

}
