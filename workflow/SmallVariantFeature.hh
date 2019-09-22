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

#include <cassert>
#include <vector>

#include "classification/SmallVariantAlignmentClassifier.hh"
#include "workflow/GraphFeature.hh"

namespace ehunter
{

class SmallVariantFeature : public GraphFeature
{
    ~SmallVariantFeature() override = default;

    explicit SmallVariantFeature(std::shared_ptr<GraphModel> modelPtr, std::vector<graphtools::NodeId> nodeIds);

    void
    process(const Read& read, const Alignments& readAligns, const Read& mate, const Alignments& mateAligns) override;

    const std::vector<ReadSummaryForSmallVariant>& readSummaries() const { return readSummaries_; }

private:
    void processRead(const Read& read, const std::list<graphtools::GraphAlignment>& alignments);

    SmallVariantAlignmentClassifier alignmentClassifier_;
    std::vector<ReadSummaryForSmallVariant> readSummaries_;
};

}
