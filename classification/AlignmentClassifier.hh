//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
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

#include <cstdint>
#include <iostream>
#include <list>
#include <set>
#include <vector>

#include "graphalign/GraphAlignment.hh"

using graphtools::GraphAlignment;

namespace ehunter
{

enum class AlignmentType
{
    kSpansRepeat,
    kFlanksRepeat,
    kInsideRepeat,
    kOutsideRepeat,
    kUnableToAlign,
    kUnprocessed
};

std::ostream& operator<<(std::ostream& os, const AlignmentType& read_class);

class RepeatAlignmentClassifier
{
public:
    RepeatAlignmentClassifier(const graphtools::Graph& graph, int32_t repeat_node_id);
    AlignmentType Classify(const GraphAlignment& alignment) const;
    GraphAlignment GetCanonicalAlignment(const std::list<GraphAlignment>& alignments) const;
    const std::set<graphtools::NodeId>& leftFlankNodeIds() const { return left_flank_node_ids_; }
    const std::set<graphtools::NodeId>& rightFlankNodeIds() const { return right_flank_node_ids_; }

    bool operator==(const RepeatAlignmentClassifier& other) const;

private:
    int32_t repeat_node_id_;
    std::set<graphtools::NodeId> left_flank_node_ids_;
    std::set<graphtools::NodeId> right_flank_node_ids_;
};

}
