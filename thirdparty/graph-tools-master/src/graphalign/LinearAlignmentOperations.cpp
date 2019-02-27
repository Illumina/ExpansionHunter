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

#include "graphalign/LinearAlignmentOperations.hh"

#include <iostream>
#include <list>
#include <sstream>
#include <stdexcept>
#include <string>

#include "graphalign/OperationOperations.hh"

using std::list;
using std::logic_error;
using std::string;

namespace graphtools
{

bool checkConsistency(const Alignment& alignment, const string& reference, const string& query)
{
    const bool does_alignment_span_whole_query = alignment.queryLength() == query.length();
    const bool is_alignment_within_reference
        = alignment.referenceStart() + alignment.referenceLength() <= reference.length();

    if (!does_alignment_span_whole_query || !is_alignment_within_reference)
    {
        return false;
    }

    uint32_t query_pos = 0;
    uint32_t ref_pos = alignment.referenceStart();

    for (const auto& operation : alignment)
    {
        const string query_piece = query.substr(query_pos, operation.queryLength());
        const string reference_piece = reference.substr(ref_pos, operation.referenceLength());

        if (!checkConsistency(operation, reference_piece, query_piece))
        {
            return false;
        }

        query_pos += operation.queryLength();
        ref_pos += operation.referenceLength();
    }

    return true;
}

list<StringPair> getSequencesForEachOperation(const Alignment& alignment, const string& reference, const string& query)
{
    list<StringPair> sequence_pieces;

    uint32_t query_pos = 0;
    uint32_t ref_pos = alignment.referenceStart();

    for (const auto& operation : alignment)
    {
        const string reference_piece = reference.substr(ref_pos, operation.referenceLength());
        const string query_piece = query.substr(query_pos, operation.queryLength());
        sequence_pieces.emplace_back(reference_piece, query_piece);

        query_pos += operation.queryLength();
        ref_pos += operation.referenceLength();
    }

    return sequence_pieces;
}

bool checkIfBookended(const Alignment& first_alignment, const Alignment& second_alignment)
{
    const size_t position_after_first_alignment_end
        = first_alignment.referenceStart() + first_alignment.referenceLength();
    const bool are_adjacent = position_after_first_alignment_end == (size_t)second_alignment.referenceStart();

    if (!are_adjacent)
    {
        return false;
    }

    const bool is_first_alignment_ends_with_softclip = first_alignment.back().type() == OperationType::kSoftclip;
    const bool is_second_alignment_starts_with_softclip = second_alignment.front().type() == OperationType::kSoftclip;

    if (is_first_alignment_ends_with_softclip || is_second_alignment_starts_with_softclip)
    {
        const bool is_first_alignment_fully_softclipped
            = first_alignment.numOperations() == 1 && is_first_alignment_ends_with_softclip;
        const bool is_second_alignment_fully_softclipped
            = second_alignment.numOperations() == 1 && is_second_alignment_starts_with_softclip;

        return is_first_alignment_fully_softclipped || is_second_alignment_fully_softclipped;
    }

    return true;
}

Alignment mergeAlignments(const Alignment& first_alignment, const Alignment& second_alignment)
{
    if (!checkIfBookended(first_alignment, second_alignment))
    {
        std::ostringstream msg;
        msg << "Alignments " << first_alignment << " and " << second_alignment << " are not bookended";
        throw std::logic_error(msg.str());
    }

    list<Operation> first_alignment_operations = first_alignment.operations();
    list<Operation> second_alignment_operations = second_alignment.operations();

    if (first_alignment_operations.back().type() == second_alignment_operations.front().type())
    {
        const Operation last_operation_of_first_alignment = first_alignment_operations.back();
        first_alignment_operations.pop_back();

        const Operation first_operations_of_second_alignment = second_alignment_operations.front();
        second_alignment_operations.pop_front();

        uint32_t merged_operation_length
            = last_operation_of_first_alignment.length() + first_operations_of_second_alignment.length();

        first_alignment_operations.emplace_back(last_operation_of_first_alignment.type(), merged_operation_length);
    }

    first_alignment_operations.splice(first_alignment_operations.end(), second_alignment_operations);

    return Alignment(first_alignment.referenceStart(), first_alignment_operations);
}

int32_t scoreAlignment(const Alignment& alignment, int32_t match_score, int32_t mismatch_score, int32_t gap_score)
{
    int32_t score = 0;

    for (const Operation& operation : alignment)
    {
        switch (operation.type())
        {
        case OperationType::kMatch:
            score += match_score * operation.referenceLength();
            break;
        case OperationType::kMismatch:
            score += mismatch_score * operation.referenceLength();
            break;
        case OperationType::kInsertionToRef:
            score += gap_score * operation.queryLength();
            break;
        case OperationType::kDeletionFromRef:
            score += gap_score * operation.referenceLength();
            break;
        default:
            break;
        }
    }

    return score;
}

string prettyPrint(const Alignment& alignment, const string& reference, const string& query)
{
    string reference_encoding, match_patten, query_encoding;

    const list<StringPair> sequence_pieces = getSequencesForEachOperation(alignment, reference, query);
    auto matched_sequences_it = sequence_pieces.begin();

    for (const auto& operation : alignment)
    {
        const string& operation_reference = matched_sequences_it->first;
        const string& operation_query = matched_sequences_it->second;

        switch (operation.type())
        {
        case OperationType::kDeletionFromRef:
            reference_encoding += operation_reference;
            query_encoding += string(operation.referenceLength(), '-');
            match_patten += string(operation.referenceLength(), ' ');
            break;

        case OperationType::kInsertionToRef:
            reference_encoding += string(operation.queryLength(), '-');
            match_patten += string(operation.queryLength(), ' ');
            query_encoding += operation_query;
            break;

        case OperationType::kMatch:
            reference_encoding += operation_reference;
            match_patten += string(operation.queryLength(), '|');
            query_encoding += operation_query;
            break;

        case OperationType::kMismatch:
        case OperationType::kMissingBases:
            reference_encoding += operation_reference;
            match_patten += string(operation.queryLength(), ' ');
            query_encoding += operation_query;
            break;

        case OperationType::kSoftclip:
            reference_encoding += string(operation.queryLength(), '-');
            match_patten += string(operation.queryLength(), ' ');
            query_encoding += operation_query;
            break;
        }

        ++matched_sequences_it;
    }

    return reference_encoding + '\n' + match_patten + '\n' + query_encoding;
}
}
