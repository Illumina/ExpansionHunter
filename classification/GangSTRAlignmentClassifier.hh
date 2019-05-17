//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Nima Mousavi <mousavi@ucsd.edu>
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

#include <cstdint>
#include <iostream>
#include <list>
#include <set>
#include <vector>

#include "graphalign/GraphAlignment.hh"

using graphtools::GraphAlignment;

namespace ehunter
{

    enum class GangSTRAlignmentType
    {
        kSpansRepeat,
        kFlanksLeft,
        kFlanksRight,
        kInsideRepeat,
        kLeftOfRepeat,
        kRightOfRepeat,
        kUnableToAlign,
        kUnprocessed
    };

    std::ostream& operator<<(std::ostream& os, const GangSTRAlignmentType& read_class);

    class GangSTRAlignmentClassifier
    {
    public:
        GangSTRAlignmentClassifier(const graphtools::Graph& graph, int32_t repeat_node_id);
        GangSTRAlignmentType Classify(const GraphAlignment& alignment) const;
        GraphAlignment GetCanonicalAlignment(const std::list<GraphAlignment>& alignments) const;
        const std::set<graphtools::NodeId>& leftFlankNodeIds() const { return left_flank_node_ids_; }
        const std::set<graphtools::NodeId>& rightFlankNodeIds() const { return right_flank_node_ids_; }

        bool operator==(const GangSTRAlignmentClassifier& other) const;

    private:
        int32_t repeat_node_id_;
        std::set<graphtools::NodeId> left_flank_node_ids_;
        std::set<graphtools::NodeId> right_flank_node_ids_;
    };

}
