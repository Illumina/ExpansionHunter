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

#pragma once

#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "thirdparty/spdlog/fmt/ostr.h"
#include "thirdparty/spdlog/spdlog.h"

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "classification/AlignmentClassifier.hh"
#include "genotyping/RepeatGenotype.hh"
#include "reads/Read.hh"
#include "region_analysis/VariantAnalyzer.hh"

namespace ehunter
{

class RepeatAnalyzer : public VariantAnalyzer
{
public:
    RepeatAnalyzer(
        std::string variantId, AlleleCount expectedAlleleCount, const graphtools::Graph& graph,
        graphtools::NodeId repeatNodeId, double haplotypeDepth, int32_t maxNumUnitsInRead)
        : VariantAnalyzer(std::move(variantId), expectedAlleleCount, graph, { repeatNodeId })
        , maxNumUnitsInRead_(maxNumUnitsInRead)
        , haplotypeDepth_(haplotypeDepth)
        , repeatUnit_(graph.nodeSeq(repeatNodeId))
        , alignmentClassifier_(graph, repeatNodeId)
    {
        verboseLogger_ = spdlog::get("verbose");
    }

    ~RepeatAnalyzer() = default;

    void processMates(
        const reads::Read& read, const graphtools::GraphAlignment& readAlignment, const reads::Read& mate,
        const graphtools::GraphAlignment& mateAlignment) override;

    std::unique_ptr<VariantFindings> analyze() const override;

private:
    graphtools::NodeId repeatNodeId() const { return nodeIds_.front(); }
    boost::optional<RepeatGenotype> genotypeCommonRepeat(const std::vector<int32_t>& candidateAlleleSizes) const;
    reads::RepeatAlignmentStats classifyReadAlignment(const graphtools::GraphAlignment& alignment);
    void summarizeAlignmentsToReadCounts(const reads::RepeatAlignmentStats& repeatAlignmentStats);
    std::vector<int32_t> generateCandidateAlleleSizes() const;

    const int32_t maxNumUnitsInRead_;
    const double haplotypeDepth_;

    const std::string repeatUnit_;
    RepeatAlignmentClassifier alignmentClassifier_;

    CountTable countsOfSpanningReads_;
    CountTable countsOfFlankingReads_;
    CountTable countsOfInrepeatReads_;

    std::shared_ptr<spdlog::logger> verboseLogger_;
};

}
