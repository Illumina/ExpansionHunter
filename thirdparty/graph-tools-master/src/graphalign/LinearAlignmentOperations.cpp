// -*- mode: c++; indent-tabs-mode: nil; -*-
//
// Copyright (c) 2018 Illumina, Inc.
// All rights reserved.

// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:

// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.

// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.

// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
        // std::cerr << "does_alignment_span_whole_query=" << does_alignment_span_whole_query << std::endl;
        // std::cerr << "is_alignment_within_reference=" << is_alignment_within_reference << std::endl;
        // std::cerr << "alignment.referenceStart()=" << alignment.referenceStart() << std::endl;
        // std::cerr << "alignment.referenceLength()=" << alignment.referenceLength() << std::endl;
        // std::cerr << "reference.length()=" << reference.length() << std::endl;
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
            // std::cerr << "operation=" << operation << std::endl;
            // std::cerr << "reference_piece=" << reference_piece << std::endl;
            // std::cerr << "query_piece=" << query_piece << std::endl;

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
