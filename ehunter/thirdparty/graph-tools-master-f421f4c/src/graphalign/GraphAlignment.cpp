//
// GraphTools library
// Copyright 2017-2019 Illumina, Inc.
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

#include <sstream>
#include <stdexcept>
#include <typeindex>

#include "graphalign/GraphAlignment.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphcore/PathOperations.hh"

using std::list;
using std::map;
using std::string;
using std::to_string;

namespace graphtools
{
void GraphAlignment::assertValidity() const
{
    for (size_t node_index = 0; node_index != path_.numNodes(); ++node_index)
    {
        const Alignment& node_alignment = alignments_[node_index];

        const bool is_start_wrong
            = path_.getStartPositionOnNodeByIndex(node_index) != (int32_t)node_alignment.referenceStart();

        const int32_t node_alignment_end = node_alignment.referenceLength() + node_alignment.referenceStart();
        const bool is_end_wrong = path_.getEndPositionOnNodeByIndex(node_index) != node_alignment_end;

        if (is_start_wrong || is_end_wrong)
        {
            std::ostringstream graph_alignment_encoding;
            graph_alignment_encoding << *this;
            throw std::logic_error(
                "Path " + path_.encode() + " is not compatible with graph alignment " + graph_alignment_encoding.str());
        }
    }
}

uint32_t GraphAlignment::queryLength() const
{
    uint32_t query_span = 0;
    for (const auto& alignment : alignments_)
    {
        query_span += alignment.queryLength();
    }
    return query_span;
}

uint32_t GraphAlignment::referenceLength() const
{
    uint32_t reference_span = 0;
    for (const auto& alignment : alignments_)
    {
        reference_span += alignment.referenceLength();
    }
    return reference_span;
}

uint32_t GraphAlignment::numMatches() const
{
    uint32_t num_matches = 0;
    for (const auto& alignment : alignments_)
    {
        num_matches += alignment.numMatched();
    }
    return num_matches;
}

bool GraphAlignment::overlapsNode(NodeId node_id) const
{
    return path_.checkOverlapWithNode(static_cast<NodeId>(node_id));
}

list<int32_t> GraphAlignment::getIndexesOfNode(NodeId node_id) const
{
    list<int32_t> indexes;
    const auto num_alignments = static_cast<int32_t>(alignments_.size());
    for (int32_t node_index = 0; node_index != num_alignments; ++node_index)
    {
        const NodeId cur_node_id = path_.getNodeIdByIndex(static_cast<size_t>(node_index));
        if (cur_node_id == node_id)
        {
            indexes.push_back(node_index);
        }
    }
    return indexes;
}

string GraphAlignment::generateCigar() const
{
    string graph_cigar;
    for (int32_t index = 0; index != (int32_t)size(); ++index)
    {
        const int32_t node_id = path_.getNodeIdByIndex(static_cast<size_t>(index));
        graph_cigar += std::to_string(node_id);
        const Alignment& alignment = alignments_[index];
        graph_cigar += "[" + alignment.generateCigar() + "]";
    }
    return graph_cigar;
}

bool GraphAlignment::operator<(const GraphAlignment& other) const
{
    if (!(path_ == other.path_))
    {
        return path_ < other.path_;
    }

    return alignments_ < other.alignments_;
}

void GraphAlignment::shrinkStart(int reference_length)
{
    int prefix_query_length = 0;

    if (reference_length >= (int)referenceLength())
    {
        std::stringstream string_stream;
        string_stream << "Cannot shrink start of " << *this << " by " << reference_length;
        throw std::logic_error(string_stream.str());
    }

    path_.shrinkStartBy(reference_length);

    size_t leftover_reference_length = reference_length;

    auto first_suffix_alignment_iter = alignments_.begin();
    while (leftover_reference_length >= first_suffix_alignment_iter->referenceLength())
    {
        leftover_reference_length -= first_suffix_alignment_iter->referenceLength();
        prefix_query_length += first_suffix_alignment_iter->queryLength();
        ++first_suffix_alignment_iter;
    }

    if (leftover_reference_length != 0)
    {
        const int split_position = first_suffix_alignment_iter->referenceStart() + leftover_reference_length;
        Alignment suffix = first_suffix_alignment_iter->splitAtReferencePosition(split_position);
        prefix_query_length += first_suffix_alignment_iter->queryLength();
        *first_suffix_alignment_iter = suffix;
    }

    Alignment softclipAlignment(first_suffix_alignment_iter->referenceStart(), to_string(prefix_query_length) + "S");
    *first_suffix_alignment_iter = mergeAlignments(softclipAlignment, *first_suffix_alignment_iter);

    alignments_.erase(alignments_.begin(), first_suffix_alignment_iter);

    assertValidity();
}

void GraphAlignment::shrinkEnd(int reference_length)
{
    int suffix_query_length = 0;

    if (reference_length >= (int)referenceLength())
    {
        std::stringstream string_stream;
        string_stream << "Cannot shrink start of " << *this << " by " << reference_length;
        throw std::logic_error(string_stream.str());
    }

    path_.shrinkEndBy(reference_length);

    size_t leftover_reference_length = reference_length;

    auto last_prefix_alignment_iter = alignments_.end() - 1;
    while (leftover_reference_length >= last_prefix_alignment_iter->referenceLength())
    {
        leftover_reference_length -= last_prefix_alignment_iter->referenceLength();
        suffix_query_length += last_prefix_alignment_iter->queryLength();
        --last_prefix_alignment_iter;
    }

    if (leftover_reference_length != 0)
    {
        const int split_position = last_prefix_alignment_iter->referenceStart()
            + last_prefix_alignment_iter->referenceLength() - leftover_reference_length;
        Alignment suffix = last_prefix_alignment_iter->splitAtReferencePosition(split_position);
        suffix_query_length += suffix.queryLength();
    }

    auto last_prefix_alignment_reference_end
        = last_prefix_alignment_iter->referenceStart() + last_prefix_alignment_iter->referenceLength();
    Alignment softclip_alignment(last_prefix_alignment_reference_end, to_string(suffix_query_length) + "S");
    *last_prefix_alignment_iter = mergeAlignments(*last_prefix_alignment_iter, softclip_alignment);

    alignments_.erase(last_prefix_alignment_iter + 1, alignments_.end());

    assertValidity();
}

std::ostream& operator<<(std::ostream& out, const GraphAlignment& graph_alignment)
{
    for (int32_t node_index = 0; node_index != (int32_t)graph_alignment.size(); ++node_index)
    {
        const int32_t node_id = graph_alignment.getNodeIdByIndex(node_index);
        out << node_id << "[" << graph_alignment[node_index] << "]";
    }
    return out;
}
}
