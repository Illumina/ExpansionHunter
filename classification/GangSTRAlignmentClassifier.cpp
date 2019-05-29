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

#include <map>
#include <string>
#include <unordered_set>

#include "classification/GangSTRAlignmentClassifier.hh"

using graphtools::NodeId;
using graphtools::Path;
using std::list;
using std::ostream;
using std::string;
using std::unordered_set;
using std::vector;

namespace ehunter
{

    ostream& operator<<(ostream& os, const GangSTRAlignmentType& read_class)
    {
        switch (read_class)
        {
            case GangSTRAlignmentType::kSpansRepeat:
                os << "kSpansRepeat";
                break;
            case GangSTRAlignmentType::kFlanksLeft:
                os << "kFlanksLeft";
                break;
            case GangSTRAlignmentType::kFlanksRight:
                os << "kFlanksRight";
                break;
            case GangSTRAlignmentType::kInsideRepeat:
                os << "kInsideRepeat";
                break;
            case GangSTRAlignmentType::kLeftOfRepeat:
                os << "kLeftOfRepeat";
                break;
            case GangSTRAlignmentType::kRightOfRepeat:
                os << "kRightOfRepeat";
                break;
            case GangSTRAlignmentType::kUnableToAlign:
                os << "kUnableToAlign";
                break;
            case GangSTRAlignmentType::kUnprocessed:
                os << "kUnprocessed";
                break;
            default:
                throw std::logic_error("Encountered unknown read class type");
        }
        return os;
    }

    GangSTRAlignmentClassifier::GangSTRAlignmentClassifier(const graphtools::Graph& graph, int32_t repeat_node_id)
            : repeat_node_id_(repeat_node_id)
    {
        left_flank_node_ids_ = graph.predecessors(repeat_node_id_);
        left_flank_node_ids_.erase(repeat_node_id_);

        right_flank_node_ids_ = graph.successors(repeat_node_id_);
        right_flank_node_ids_.erase(repeat_node_id_);
    }

    GraphAlignment GangSTRAlignmentClassifier::GetCanonicalAlignment(const list<GraphAlignment>& Alignments) const
    {
        const GraphAlignment* canonical_alignment_ptr = nullptr;
        for (const GraphAlignment& alignment : Alignments)
        {
            GangSTRAlignmentType alignment_type = Classify(alignment);
            if (!canonical_alignment_ptr)
            {
                canonical_alignment_ptr = &alignment;
            }

            if (alignment_type == GangSTRAlignmentType::kInsideRepeat)
            {
                return alignment;
            }
            // TODO Nima: Make sure this achieves what is intended
            else if (alignment_type == GangSTRAlignmentType::kFlanksRight ||
                    alignment_type == GangSTRAlignmentType::kFlanksLeft)
            {
                canonical_alignment_ptr = &alignment;
            }
        }
        return *canonical_alignment_ptr;
    }

    GangSTRAlignmentType GangSTRAlignmentClassifier::Classify(const GraphAlignment& alignment) const
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
//        const bool overlaps_either_flank = overlaps_left_flank || overlaps_right_flank;

        if (overlaps_both_flanks)
        {
            return GangSTRAlignmentType::kSpansRepeat;
        }

        const bool overlaps_repeat = alignment.overlapsNode(repeat_node_id_);
        if (overlaps_left_flank)
        {
            if (overlaps_repeat) {
                return GangSTRAlignmentType::kFlanksLeft;
            }
            else {
                return GangSTRAlignmentType::kLeftOfRepeat;
            }
        }
        if (overlaps_right_flank)
        {
            if (overlaps_repeat){
                return GangSTRAlignmentType::kFlanksRight;
            }
            else {
                return GangSTRAlignmentType::kRightOfRepeat;
            }
        }
        if (overlaps_repeat)
        {
            return GangSTRAlignmentType::kInsideRepeat;
        }

        // TODO Nima Is this class appropriate? This is a default case that shouldn't happen.
        return GangSTRAlignmentType::kUnableToAlign;
    }

    bool GangSTRAlignmentClassifier::operator==(const GangSTRAlignmentClassifier& other) const
    {
        return repeat_node_id_ == other.repeat_node_id_ && left_flank_node_ids_ == other.left_flank_node_ids_
               && right_flank_node_ids_ == other.right_flank_node_ids_;
    }

}
