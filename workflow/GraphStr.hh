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

#include "graphcore/Graph.hh"

#include "classification/AlignmentSummary.hh"
#include "classification/StrAlignmentClassifier.hh"
#include "strs/StrAlignmentStats.hh"
#include "workflow/GraphFeature.hh"

namespace ehunter
{

class GraphModel;

class GraphStr : public GraphFeature
{
public:
    GraphStr(std::shared_ptr<GraphModel> model, graphtools::NodeId motifNodeId);
    ~GraphStr() override = default;
    std::shared_ptr<RegionModel> model() override;

    void summarize(
        const std::string& read, const Alignments& readAligns, const std::string& mate,
        const Alignments& mateAligns) override;

    const std::string& motif() const;
    const std::vector<ReadSummaryForStr>& readSummaries() const { return readSummaries_; }
    StrAlignmentStats alignmentStats() const { return statsCalculator_.getStats(); }

private:
    std::shared_ptr<GraphModel> model_;
    graphtools::NodeId motifNode_;

    StrAlignmentClassifier alignmentClassifier_;
    StrAlignmentStatsCalculator statsCalculator_;
    std::vector<ReadSummaryForStr> readSummaries_;
};

}
