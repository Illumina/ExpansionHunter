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
        graphtools::NodeId repeatNodeId)
        : VariantAnalyzer(std::move(variantId), expectedAlleleCount, graph, { repeatNodeId })
        , repeatUnit_(graph.nodeSeq(repeatNodeId))
        , alignmentClassifier_(graph, repeatNodeId)
        , console_(spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console"))
    {
    }

    ~RepeatAnalyzer() = default;

    void processMates(
        const Read& read, const graphtools::GraphAlignment& readAlignment, const Read& mate,
        const graphtools::GraphAlignment& mateAlignment) override;

    std::unique_ptr<VariantFindings> analyze(const SampleParameters& params) const override;

private:
    graphtools::NodeId repeatNodeId() const { return nodeIds_.front(); }
    RepeatAlignmentStats classifyReadAlignment(const graphtools::GraphAlignment& alignment);
    void summarizeAlignmentsToReadCounts(const RepeatAlignmentStats& repeatAlignmentStats);
    std::vector<int32_t> generateCandidateAlleleSizes() const;

    const std::string repeatUnit_;
    RepeatAlignmentClassifier alignmentClassifier_;

    CountTable countsOfSpanningReads_;
    CountTable countsOfFlankingReads_;
    CountTable countsOfInrepeatReads_;

    std::shared_ptr<spdlog::logger> console_;
};

}
