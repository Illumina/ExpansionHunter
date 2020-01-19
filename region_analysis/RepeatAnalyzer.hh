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
#include "filtering/GraphVariantAlignmentStats.hh"
#include "genotyping/RepeatGenotype.hh"
#include "reads/Read.hh"
#include "region_analysis/VariantAnalyzer.hh"

namespace ehunter
{

class RepeatAnalyzer : public VariantAnalyzer
{
public:
    RepeatAnalyzer(
        std::string variantId, const graphtools::Graph& graph, graphtools::NodeId repeatNodeId,
        GenotyperParameters genotyperParams)
        : VariantAnalyzer(std::move(variantId), graph, { repeatNodeId }, genotyperParams)
        , repeatUnit_(graph.nodeSeq(repeatNodeId))
        , alignmentClassifier_(graph, repeatNodeId)
        , alignmentStatsCalculator_({ repeatNodeId })
        , countOfInrepeatReadPairs_(0)
        , console_(spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console"))
    {
    }

    ~RepeatAnalyzer() = default;

    const std::string& repeatUnit() const { return repeatUnit_; }
    void addInrepeatReadPair() { countOfInrepeatReadPairs_++; }

    void processMates(
        const Read& read, const graphtools::GraphAlignment& readAlignment, const Read& mate,
        const graphtools::GraphAlignment& mateAlignment) override;

    std::unique_ptr<VariantFindings> analyze(const LocusStats& stats) const override;

private:
    graphtools::NodeId repeatNodeId() const { return nodeIds_.front(); }
    RepeatAlignmentStats classifyReadAlignment(const graphtools::GraphAlignment& alignment);
    void summarizeAlignmentsToReadCounts(const RepeatAlignmentStats& repeatAlignmentStats);

    GenotyperParameters genotyperParams_;
    const std::string repeatUnit_;
    RepeatAlignmentClassifier alignmentClassifier_;
    GraphVariantAlignmentStatsCalculator alignmentStatsCalculator_;

    CountTable countsOfSpanningReads_;
    CountTable countsOfFlankingReads_;
    CountTable countsOfInrepeatReads_;
    int countOfInrepeatReadPairs_;

    std::shared_ptr<spdlog::logger> console_;
};

}
