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

#include "workflow/GraphStr.hh"
#include "workflow/GraphModel.hh"

using std::shared_ptr;
using std::static_pointer_cast;

namespace ehunter
{

GraphStr::GraphStr(std::shared_ptr<GraphModel> model, graphtools::NodeId nodeId)
    : GraphVariant(std::move(model), { nodeId })
    , alignmentClassifier_(model_->graph(), nodeId)
{
}

void GraphStr::process(const Read& read, const Alignments& readAligns, const Read& mate, const Alignments& mateAligns)
{
    ReadSummaryForStr strRead = alignmentClassifier_.classifyRead(read.sequence(), readAligns);
    if (strRead.hasAlignments())
    {
        readSummaries_.push_back(std::move(strRead));
    }

    ReadSummaryForStr strMate = alignmentClassifier_.classifyRead(mate.sequence(), mateAligns);
    if (strMate.hasAlignments())
    {
        readSummaries_.push_back(std::move(strMate));
    }
}

graphtools::NodeId GraphStr::motifNodeId() const { return nodeIds_.front(); }

const std::string& GraphStr::motif() const { return model_->graph().nodeSeq(motifNodeId()); }

}
