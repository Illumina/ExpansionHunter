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

#include <string>
#include <vector>

#include "region/GraphRegion.hh"

#include "classification/AlignmentSummary.hh"
#include "classification/StrAlignmentClassifier.hh"

namespace ehunter
{

class StrFeature : public GraphFeature
{
public:
    explicit StrFeature(GraphRegion::SPtr graphModelPtr, graphtools::NodeId nodeId)
        : GraphFeature(graphModelPtr, { nodeId })
        , alignmentClassifier_(graphModelPtr->graph(), nodeId)
    {
    }

    void
    process(const Read& read, const Alignments& readAligns, const Read& mate, const Alignments& mateAligns) override;

    const std::string& motif() const;
    const std::vector<ReadSummaryForStr>& readSummaries() const { return readSummaries_; }

private:
    StrAlignmentClassifier alignmentClassifier_;
    std::vector<ReadSummaryForStr> readSummaries_;
};

}
