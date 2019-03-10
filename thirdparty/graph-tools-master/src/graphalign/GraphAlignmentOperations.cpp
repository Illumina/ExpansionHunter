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

#include "graphalign/GraphAlignmentOperations.hh"

#include <stdexcept>

#include <boost/range/adaptor/reversed.hpp>

#include "graphalign/LinearAlignmentOperations.hh"

using std::list;
using std::string;
using std::vector;

namespace graphtools
{

bool checkConsistency(const GraphAlignment& graph_alignment, const string& query)
{
    NodeId node_index = 0;
    uint32_t query_pos = 0;
    const Graph* graph_ptr = graph_alignment.path().graphRawPtr();
    for (const auto& linear_alignment : graph_alignment)
    {
        if (query.length() < query_pos + linear_alignment.queryLength())
        {
            return false;
        }
        const string query_piece = query.substr(query_pos, linear_alignment.queryLength());
        query_pos += linear_alignment.queryLength();

        NodeId node_id = graph_alignment.getNodeIdByIndex(node_index);
        const string& node_seq = graph_ptr->nodeSeq(node_id);
        if (!checkConsistency(linear_alignment, node_seq, query_piece))
        {
            return false;
        }

        const bool query_or_reference_length_is_zero
            = linear_alignment.referenceLength() == 0 || linear_alignment.queryLength() == 0;

        if (graph_alignment.path().numNodes() != 1 && query_or_reference_length_is_zero)
        {
            return false;
        }

        ++node_index;
    }

    return true;
}

static bool startsWithMatch(const Alignment& alignment)
{
    for (const auto& operation : alignment)
    {
        if (operation.type() != OperationType::kSoftclip)
        {
            return operation.type() == OperationType::kMatch;
        }
    }

    return false;
}

static bool endsWithMatch(const Alignment& alignment)
{
    for (const auto& operation : boost::adaptors::reverse(alignment))
    {
        if (operation.type() != OperationType::kSoftclip)
        {
            return operation.type() == OperationType::kMatch;
        }
    }

    return false;
}

bool isLocalAlignment(const GraphAlignment& graph_alignment)
{
    return startsWithMatch(graph_alignment.front()) && endsWithMatch(graph_alignment.back());
}

static vector<string> splitGraphCigar(const string& graph_cigar)
{
    vector<string> node_cigars;
    string node_cigar;
    for (size_t index = 0; index != graph_cigar.length(); ++index)
    {
        node_cigar += graph_cigar[index];
        if (node_cigar.back() == ']')
        {
            node_cigars.push_back(node_cigar);
            node_cigar.clear();
        }
    }

    return node_cigars;
}

GraphAlignment decodeGraphAlignment(int32_t first_node_start, const string& graph_cigar, const Graph* graph_ptr)
{
    vector<NodeId> node_ids;
    vector<Alignment> alignments;
    vector<string> node_cigars = splitGraphCigar(graph_cigar);
    for (const string& node_cigar : node_cigars)
    {
        int32_t ref_pos = alignments.empty() ? first_node_start : 0;

        string cigar;
        NodeId node_id;
        splitNodeCigar(node_cigar, cigar, node_id);
        node_ids.push_back(node_id);

        Alignment alignment(ref_pos, cigar);
        alignments.push_back(alignment);
    }

    // Convert to 0-based coordinates
    int32_t last_node_end = alignments.back().referenceStart() + alignments.back().referenceLength();
    Path path(graph_ptr, first_node_start, node_ids, last_node_end);
    return GraphAlignment(path, alignments);
}

void splitNodeCigar(const string& node_cigar, string& cigar, NodeId& node_id)
{
    node_id = static_cast<NodeId>(-1);
    string nodeid_encoding;
    for (size_t index = 0; index != node_cigar.length(); ++index)
    {
        if (node_cigar[index] == '[')
        {
            node_id = static_cast<NodeId>(std::stoull(nodeid_encoding));
            cigar = node_cigar.substr(index + 1);
            cigar.pop_back();
            return;
        }
        if (isdigit(node_cigar[index]) == 0)
        {
            throw std::logic_error(node_cigar + " is a malformed node CIGAR");
        }
        nodeid_encoding += node_cigar[index];
    }
}

// Implementation note: Unless stated otherwise, all calculations are performed in path coordinates.
GraphAlignment projectAlignmentOntoGraph(Alignment linear_alignment, Path path)
{
    vector<Alignment> alignments;

    path.shrinkBy(
        linear_alignment.referenceStart(),
        path.length() - linear_alignment.referenceStart() - linear_alignment.referenceLength());
    linear_alignment.setReferenceStart(0);

    for (size_t node_index = 0; node_index != path.numNodes(); ++node_index)
    {
        const size_t last_position_of_path_on_this_node = path.getNodeOverlapLengthByIndex(node_index);

        const size_t linear_alignment_last_position
            = linear_alignment.referenceStart() + linear_alignment.referenceLength();

        if (linear_alignment_last_position <= last_position_of_path_on_this_node)
        {
            alignments.push_back(linear_alignment);
            break;
        }
        else
        {
            Alignment suffix = linear_alignment.splitAtReferencePosition(last_position_of_path_on_this_node);
            alignments.push_back(linear_alignment);
            linear_alignment = suffix;
            linear_alignment.setReferenceStart(0);
        }
    }

    Alignment& first_alignment = alignments.front();
    first_alignment.setReferenceStart(path.startPosition());

    return GraphAlignment(path, alignments);
}

list<string> getQuerySequencesForEachNode(const GraphAlignment& graph_alignment, const string& query)
{
    list<string> sequence_pieces;

    uint32_t query_pos = 0;
    for (const auto& linear_alignment : graph_alignment)
    {
        const string query_piece = query.substr(query_pos, linear_alignment.queryLength());
        sequence_pieces.push_back(query_piece);
        query_pos += linear_alignment.queryLength();
    }

    return sequence_pieces;
}

static string joinLinearAlignmentEncodings(const vector<string>& linear_alignment_encodings)
{
    string graph_query_encoding, graph_match_pattern_encoding, graph_reference_encoding;

    for (const string& linear_alignment_encoding : linear_alignment_encodings)
    {
        if (!graph_query_encoding.empty())
        {
            graph_query_encoding += ':';
            graph_match_pattern_encoding += ':';
            graph_reference_encoding += ':';
        }

        const vector<string> query_pattern_reference_encodings
            = splitStringByDelimiter(linear_alignment_encoding, '\n');

        graph_query_encoding += query_pattern_reference_encodings[0];
        graph_match_pattern_encoding += query_pattern_reference_encodings[1];
        graph_reference_encoding += query_pattern_reference_encodings[2];
    }

    return graph_query_encoding + '\n' + graph_match_pattern_encoding + '\n' + graph_reference_encoding;
}

string prettyPrint(const GraphAlignment& graph_alignment, const string& query)
{
    const list<string> node_queries = getQuerySequencesForEachNode(graph_alignment, query);

    const Graph* graph_ptr = graph_alignment.path().graphRawPtr();

    vector<string> linear_alignment_encodings;
    size_t node_index = 0;
    for (const string& node_query : node_queries)
    {
        NodeId node_id = graph_alignment.getNodeIdByIndex(node_index);
        const string& node_seq = graph_ptr->nodeSeq(node_id);

        const auto& linear_alignment = graph_alignment[node_index];
        linear_alignment_encodings.push_back(prettyPrint(linear_alignment, node_seq, node_query));

        ++node_index;
    }

    return joinLinearAlignmentEncodings(linear_alignment_encodings);
}
}
