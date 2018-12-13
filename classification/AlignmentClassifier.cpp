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

#include <map>
#include <string>
#include <unordered_set>

#include "classification/AlignmentClassifier.hh"

using graphtools::NodeId;
using graphtools::Path;
using std::list;
using std::ostream;
using std::string;
using std::unordered_set;
using std::vector;

namespace ehunter
{

ostream& operator<<(ostream& os, const AlignmentType& read_class)
{
    switch (read_class)
    {
    case AlignmentType::kSpansRepeat:
        os << "kSpansRepeat";
        break;
    case AlignmentType::kFlanksRepeat:
        os << "kFlanksRepeat";
        break;
    case AlignmentType::kInsideRepeat:
        os << "kInsideRepeat";
        break;
    case AlignmentType::kOutsideRepeat:
        os << "kOutsideRepeat";
        break;
    case AlignmentType::kUnableToAlign:
        os << "kUnableToAlign";
        break;
    case AlignmentType::kUnprocessed:
        os << "kUnprocessed";
        break;
    default:
        throw std::logic_error("Encountered unknown read class type");
    }
    return os;
}

RepeatAlignmentClassifier::RepeatAlignmentClassifier(const graphtools::Graph& graph, int32_t repeat_node_id)
    : repeat_node_id_(repeat_node_id)
{
    left_flank_node_ids_ = graph.predecessors(repeat_node_id_);
    left_flank_node_ids_.erase(repeat_node_id_);

    right_flank_node_ids_ = graph.successors(repeat_node_id_);
    right_flank_node_ids_.erase(repeat_node_id_);
}

GraphAlignment RepeatAlignmentClassifier::GetCanonicalAlignment(const list<GraphAlignment>& Alignments) const
{
    const GraphAlignment* canonical_alignment_ptr = nullptr;
    for (const GraphAlignment& alignment : Alignments)
    {
        AlignmentType alignment_type = Classify(alignment);
        if (!canonical_alignment_ptr)
        {
            canonical_alignment_ptr = &alignment;
        }

        if (alignment_type == AlignmentType::kInsideRepeat)
        {
            return alignment;
        }
        else if (alignment_type == AlignmentType::kFlanksRepeat)
        {
            canonical_alignment_ptr = &alignment;
        }
    }
    return *canonical_alignment_ptr;
}

AlignmentType RepeatAlignmentClassifier::Classify(const GraphAlignment& alignment) const
{
    bool overlaps_left_flank = false;
    bool overlaps_right_flank = false;

    for (auto node_id : alignment.path().nodeIds())
    {
        if (left_flank_node_ids_.find(node_id) != left_flank_node_ids_.end())
        {
            overlaps_left_flank = true;
        }

        if (right_flank_node_ids_.find(node_id) != right_flank_node_ids_.end())
        {
            overlaps_right_flank = true;
        }
    }

    const bool overlaps_both_flanks = overlaps_left_flank && overlaps_right_flank;
    const bool overlaps_either_flank = overlaps_left_flank || overlaps_right_flank;

    if (overlaps_both_flanks)
    {
        return AlignmentType::kSpansRepeat;
    }

    const bool overlaps_repeat = alignment.overlapsNode(repeat_node_id_);
    if (overlaps_either_flank && overlaps_repeat)
    {
        return AlignmentType::kFlanksRepeat;
    }

    if (overlaps_repeat)
    {
        return AlignmentType::kInsideRepeat;
    }

    return AlignmentType::kOutsideRepeat;
}

bool RepeatAlignmentClassifier::operator==(const RepeatAlignmentClassifier& other) const
{
    return repeat_node_id_ == other.repeat_node_id_ && left_flank_node_ids_ == other.left_flank_node_ids_
        && right_flank_node_ids_ == other.right_flank_node_ids_;
}

}
