//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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

#include "locus/LocusAnalyzer.hh"

#include <boost/smart_ptr/make_unique.hpp>

#include "locus/LocusAligner.hh"
#include "locus/RFC1MotifAnalysis.hh"
#include "locus/RepeatAnalyzer.hh"
#include "locus/SmallVariantAnalyzer.hh"

using boost::make_unique;
using boost::optional;
using graphtools::AlignmentWriter;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using std::string;
using std::vector;

namespace ehunter
{
namespace locus
{

LocusAnalyzer::LocusAnalyzer(LocusSpecification locusSpec, const HeuristicParameters& params, AlignWriterPtr writer)
    : locusSpec_(std::move(locusSpec))
    , alignmentBuffer_(locusSpec_.useRFC1MotifAnalysis() ? std::make_shared<locus::AlignmentBuffer>() : nullptr)
    , aligner_(locusSpec_.locusId(), &locusSpec_.regionGraph(), params, std::move(writer), alignmentBuffer_)
    , statsCalc_(locusSpec_.typeOfChromLocusLocatedOn(), locusSpec_.regionGraph())
{
    for (const auto& variantSpec : locusSpec_.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kRepeat)
        {
            const auto& graph = locusSpec_.regionGraph();
            const int repeatNodeId = static_cast<int>(variantSpec.nodes().front());
            const auto& motif = graph.nodeSeq(repeatNodeId);

            if (variantSpec.classification().subtype == VariantSubtype::kRareRepeat)
            {
                if (irrPairFinder())
                {
                    const string message
                        = "Region " + locusSpec_.locusId() + " must not have more than one rare repeat";
                    throw std::logic_error(message);
                }
                addIrrPairFinder(motif);
            }

            addRepeatAnalyzer(variantSpec.id(), repeatNodeId);
        }
        else if (variantSpec.classification().type == VariantType::kSmallVariant)
        {
            addSmallVariantAnalyzer(
                variantSpec.id(), variantSpec.classification().subtype, variantSpec.nodes(),
                variantSpec.optionalRefNode());
        }
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw std::logic_error("Missing logic to create an analyzer for " + encoding.str());
        }
    }
}

void LocusAnalyzer::processMates(
    Read& read, Read* mate, RegionType regionType, graphtools::AlignerSelector& alignerSelector)
{
    if (regionType == RegionType::kTarget)
    {
        processOntargetMates(read, mate, alignerSelector);
    }
    else if (mate)
    {
        processOfftargetMates(read, *mate);
    }
}

void LocusAnalyzer::processOntargetMates(Read& read, Read* mate, graphtools::AlignerSelector& alignerSelector)
{
    auto alignedPair = aligner_.align(read, mate, alignerSelector);

    const bool neitherMateAligned = !alignedPair.first && !alignedPair.second;
    const bool bothMatesAligned = alignedPair.first && alignedPair.second;

    if (irrPairFinder_ && neitherMateAligned && mate)
    {
        processOfftargetMates(read, *mate);
        return;
    }

    if (bothMatesAligned)
    {
        statsCalc_.inspect(*alignedPair.first, *alignedPair.second);
        runVariantAnalysis(read, *alignedPair.first, *mate, *alignedPair.second);
    }
    else
    {
        if (alignedPair.first)
        {
            statsCalc_.inspectRead(*alignedPair.first);
        }
        if (alignedPair.second)
        {
            statsCalc_.inspectRead(*alignedPair.second);
        }
    }
}

void LocusAnalyzer::processOfftargetMates(const Read& read, const Read& mate)
{
    if (!irrPairFinder_)
    {
        const string message = "Locus " + locusSpec_.locusId() + " is not supposed to have offtarget read pairs";
        throw std::logic_error(message);
    }

    if (irrPairFinder_->check(read.sequence(), mate.sequence()))
    {
        int numAnalyzersFound = 0;
        for (auto& variantAnalyzer : variantAnalyzers_)
        {
            auto repeatAnalyzer = dynamic_cast<RepeatAnalyzer*>(variantAnalyzer.get());
            if (repeatAnalyzer != nullptr && repeatAnalyzer->repeatUnit() == irrPairFinder_->targetMotif())
            {
                numAnalyzersFound++;
                repeatAnalyzer->addInrepeatReadPair();
            }
        }

        if (numAnalyzersFound != 1)
        {
            const string message = "Locus " + locusSpec_.locusId() + " must have exactly one rare motif";
            throw std::logic_error(message);
        }
    }
}

LocusFindings LocusAnalyzer::analyze(Sex sampleSex, boost::optional<double> genomeWideDepth)
{
    LocusFindings locusFindings(statsCalc_.estimate(sampleSex));
    if (genomeWideDepth && locusSpec_.requiresGenomeWideDepth())
    {
        locusFindings.stats.setDepth(*genomeWideDepth);
    }

    for (auto& variantAnalyzer : variantAnalyzers_)
    {
        std::unique_ptr<VariantFindings> variantFindingsPtr = variantAnalyzer->analyze(locusFindings.stats);
        const string& variantId = variantAnalyzer->variantId();
        locusFindings.findingsForEachVariant.emplace(variantId, std::move(variantFindingsPtr));
    }

    // Run RFC1 caller if required for this locus:
    if (locusSpec().useRFC1MotifAnalysis())
    {
        assert(alignmentBuffer_);
        runRFC1MotifAnalysis(*alignmentBuffer_, locusFindings);
    }

    return locusFindings;
}

void LocusAnalyzer::addIrrPairFinder(std::string motif) { irrPairFinder_ = IrrPairFinder(std::move(motif)); }

void LocusAnalyzer::addRepeatAnalyzer(std::string variantId, graphtools::NodeId nodeId)
{
    variantAnalyzers_.emplace_back(make_unique<RepeatAnalyzer>(
        std::move(variantId), locusSpec_.regionGraph(), nodeId, locusSpec_.genotyperParameters()));
}

void LocusAnalyzer::addSmallVariantAnalyzer(
    string variantId, VariantSubtype subtype, vector<NodeId> nodes, optional<NodeId> refNode)
{
    variantAnalyzers_.emplace_back(make_unique<SmallVariantAnalyzer>(
        std::move(variantId), subtype, locusSpec_.regionGraph(), std::move(nodes), refNode,
        locusSpec_.genotyperParameters()));
}

void LocusAnalyzer::runVariantAnalysis(
    const Read& read, const LocusAnalyzer::Align& readAlign, const Read& mate, const LocusAnalyzer::Align& mateAlign)
{
    for (auto& analyzer : variantAnalyzers_)
    {
        analyzer->processMates(read, readAlign, mate, mateAlign);
    }
}

}
}
