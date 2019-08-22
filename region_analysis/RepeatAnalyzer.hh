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

#include <list>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "thirdparty/spdlog/fmt/ostr.h"
#include "thirdparty/spdlog/spdlog.h"

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "classification/AlignmentSummary.hh"
#include "classification/StrAlignmentClassifier.hh"
#include "genotyping/RepeatGenotype.hh"
#include "reads/Read.hh"
#include "region_analysis/VariantAnalyzer.hh"

namespace ehunter
{

class RepeatAnalyzer : public VariantAnalyzer
{
public:
    RepeatAnalyzer(std::string variantId, const graphtools::Graph& graph, graphtools::NodeId repeatNodeId)
        : VariantAnalyzer(std::move(variantId), graph, { repeatNodeId })
        , repeatUnit_(graph.nodeSeq(repeatNodeId))
        , alignmentClassifier_(graph, repeatNodeId)
        , countOfInrepeatReadPairs_(0)
        , console_(spdlog::get("console") ? spdlog::get("console") : spdlog::stderr_color_mt("console"))
    {
    }

    ~RepeatAnalyzer() = default;

    const std::string& repeatUnit() const { return repeatUnit_; }
    void addInrepeatReadPair() { countOfInrepeatReadPairs_++; }

    void processMates(
        const Read& read, const std::list<GraphAlignment>& readAlignments, const Read& mate,
        const std::list<GraphAlignment>& mateAlignments) override;

    std::unique_ptr<VariantFindings> analyze(const LocusStats& stats) const override;

private:
    graphtools::NodeId repeatNodeId() const { return nodeIds_.front(); }

    // std::set<StrAlignment> summarize(const Read& read, const std::list<GraphAlignment>& alignments) const;

    void summarizeAlignmentsToReadCounts(const StrAlignment& strAlignment);

    std::string repeatUnit_;
    StrAlignmentClassifier alignmentClassifier_;

    CountTable countsOfSpanningReads_;
    CountTable countsOfFlankingReads_;
    CountTable countsOfInrepeatReads_;
    int countOfInrepeatReadPairs_;

    std::shared_ptr<spdlog::logger> console_;
};

}
