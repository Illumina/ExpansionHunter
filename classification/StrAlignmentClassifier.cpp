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

#include <map>
#include <set>
#include <string>
#include <unordered_set>

#include "alignment/AlignmentFilters.hh"
#include "alignment/OperationsOnAlignments.hh"
#include "classification/StrAlignmentClassifier.hh"
#include "stats/WeightedPurityCalculator.hh"

using boost::optional;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::Path;
using std::list;
using std::ostream;
using std::set;
using std::string;
using std::unordered_set;
using std::vector;

namespace ehunter
{

bool StrAlignmentClassifier::checkQuality(
    const string& read, const GraphAlignment& alignment, const StrAlignment& strAlignment) const
{

    if (!checkIfPassesAlignmentFilters(alignment))
    {
        return false;
    }

    const bool doesReadAlignWellOverLeftFlank = checkIfUpstreamAlignmentIsGood(repeatNodeId_, alignment);
    const bool doesReadAlignWellOverRightFlank = checkIfDownstreamAlignmentIsGood(repeatNodeId_, alignment);

    if (strAlignment.type() == StrAlignment::Type::kFlanking)
    {
        if (!doesReadAlignWellOverLeftFlank && !doesReadAlignWellOverRightFlank)
        {
            return false;
        }
    }

    if (strAlignment.type() == StrAlignment::Type::kSpanning)
    {
        if (!doesReadAlignWellOverLeftFlank || !doesReadAlignWellOverRightFlank)
        {
            return false;
        }
    }

    if (strAlignment.type() == StrAlignment::Type::kInrepeat)
    {
        const string& repeatUnit = alignment.path().graphRawPtr()->nodeSeq(repeatNodeId_);
        WeightedPurityCalculator wpCalculator(repeatUnit);
        const double kScoreCutoff = 0.8;
        return wpCalculator.score(read) >= kScoreCutoff;
    }

    return true;
}

StrAlignmentClassifier::StrAlignmentClassifier(const graphtools::Graph& graph, int repeatNodeId)
    : repeatNodeId_(repeatNodeId)
{
    leftFlankNodeIds_ = graph.predecessors(repeatNodeId_);
    leftFlankNodeIds_.erase(repeatNodeId_);

    rightFlankNodeIds_ = graph.successors(repeatNodeId_);
    rightFlankNodeIds_.erase(repeatNodeId_);
}

ReadSummaryForStr StrAlignmentClassifier::classifyRead(const string& read, const list<GraphAlignment>& alignments) const
{
    ReadSummaryForStr strRead(read.length());

    for (const auto& alignment : alignments)
    {
        optional<StrAlignment> optionalSummary = classify(alignment);
        if (optionalSummary && checkQuality(read, alignment, *optionalSummary))
        {
            strRead.addAlignment(*optionalSummary);
        }
    }

    return strRead;
}

optional<StrAlignment> StrAlignmentClassifier::classify(const GraphAlignment& alignment) const
{
    bool overlapsLeftFlank = false;
    bool overlapsRightFlank = false;

    for (auto node_id : alignment.path().nodeIds())
    {
        if (leftFlankNodeIds_.find(node_id) != leftFlankNodeIds_.end())
        {
            overlapsLeftFlank = true;
        }

        if (rightFlankNodeIds_.find(node_id) != rightFlankNodeIds_.end())
        {
            overlapsRightFlank = true;
        }
    }

    const bool overlapsBothFlanks = overlapsLeftFlank && overlapsRightFlank;
    const bool overlapsEitherFlank = overlapsLeftFlank || overlapsRightFlank;

    const int score = scoreAlignment(alignment);
    const int numRepeatUnitsOverlapped = countFullOverlaps(repeatNodeId_, alignment);

    if (overlapsBothFlanks)
    {
        return StrAlignment(numRepeatUnitsOverlapped, StrAlignment::Type::kSpanning, score,0);
    }

    const bool overlaps_repeat = alignment.overlapsNode(repeatNodeId_);
    if (overlapsEitherFlank && overlaps_repeat)
    {
        return StrAlignment(numRepeatUnitsOverlapped, StrAlignment::Type::kFlanking, score,0);
    }

    if (overlaps_repeat)
    {
        return StrAlignment(numRepeatUnitsOverlapped, StrAlignment::Type::kInrepeat, score,0);
    }

    return boost::none;
}

bool StrAlignmentClassifier::operator==(const StrAlignmentClassifier& other) const
{
    return repeatNodeId_ == other.repeatNodeId_ && leftFlankNodeIds_ == other.leftFlankNodeIds_
        && rightFlankNodeIds_ == other.rightFlankNodeIds_;
}

}
