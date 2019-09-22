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

#include "workflow/SmallVariantFeature.hh"

namespace ehunter
{

void SmallVariantFeature::process(
    const Read& read, const Alignments& readAligns, const Read& mate, const Alignments& mateAligns)
{
    processRead(read, readAligns);
    processRead(mate, mateAligns);
}

void SmallVariantFeature::processRead(const Read& read, const std::list<graphtools::GraphAlignment>& alignments)
{
    ReadSummaryForSmallVariant smallVariantRead = alignmentClassifier_.classifyRead(read.sequence(), alignments);

    if (smallVariantRead.numAlignments() > 0)
    {
        readSummaries_.push_back(std::move(smallVariantRead));
    }
}

SmallVariantFeature::SmallVariantFeature(std::shared_ptr<GraphModel> modelPtr, std::vector<graphtools::NodeId> nodeIds)
    : GraphFeature(modelPtr, std::move(nodeIds))
    , alignmentClassifier_(nodeIds_)

{
}

}
