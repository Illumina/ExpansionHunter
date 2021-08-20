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

#include "graphalign/LinearAlignment.hh"

#include <algorithm>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "graphalign/OperationOperations.hh"
#include "graphutils/BaseMatching.hh"

using std::list;
using std::logic_error;
using std::map;
using std::string;
using std::to_string;
using std::vector;

namespace graphtools
{

Alignment::Alignment(uint32_t reference_start, const string& cigar)
    : reference_start_(reference_start)
{
    decodeCigar(cigar);
    updateCounts();
}

void Alignment::updateCounts()
{
    clipped_ = 0;
    matched_ = 0;
    mismatched_ = 0;
    missing_ = 0;
    inserted_ = 0;
    deleted_ = 0;
    for (const Operation& operation : operations_)
    {
        switch (operation.type())
        {
        case OperationType::kSoftclip:
            clipped_ += operation.length();
            break;
        case OperationType::kMatch:
            matched_ += operation.length();
            break;
        case OperationType::kMismatch:
            mismatched_ += operation.length();
            break;
        case OperationType::kMissingBases:
            missing_ += operation.length();
            break;
        case OperationType::kInsertionToRef:
            inserted_ += operation.length();
            break;
        case OperationType::kDeletionFromRef:
            deleted_ += operation.length();
            break;
        }
    }
}

void Alignment::decodeCigar(const string& cigar)
{
    string length_encoding;
    for (char c : cigar)
    {
        if (isalpha(c) != 0)
        {
            uint32_t operation_length = std::stoi(length_encoding);
            OperationType operation_type = decodeOperationType(c);
            operations_.emplace_back(operation_type, operation_length);
            length_encoding.clear();
        }
        else
        {
            if (isdigit(c) == 0)
            {
                throw logic_error(cigar + " is malformed CIGAR string");
            }
            length_encoding += c;
        }
    }
}

uint32_t Alignment::queryLength() const
{
    int32_t query_span = 0;
    for (const auto& operation : operations_)
    {
        query_span += operation.queryLength();
    }
    return query_span;
}

uint32_t Alignment::referenceLength() const
{
    int32_t reference_span = 0;
    for (const auto& operation : operations_)
    {
        reference_span += operation.referenceLength();
    }
    return reference_span;
}

string Alignment::generateCigar() const
{
    string cigar_string;
    for (const auto& operation : operations_)
    {
        cigar_string += operation.generateCigar();
    }
    return cigar_string;
}

Alignment Alignment::splitAtReferencePosition(size_t reference_position)
{
    const size_t end_of_reference_positions = referenceStart() + referenceLength();
    if (reference_position == 0 || end_of_reference_positions <= reference_position)
    {
        std::ostringstream os;
        os << *this;
        throw logic_error("Cannot split " + os.str() + " at reference position " + to_string(reference_position));
    }

    size_t first_unused_position = reference_start_;

    list<Operation>::const_iterator operation_it = operations_.begin();

    while (operation_it != operations_.end())
    {
        const size_t first_unused_position_after_applying_operation
            = first_unused_position + operation_it->referenceLength();
        if (first_unused_position_after_applying_operation <= reference_position)
        {
            ++operation_it;
            first_unused_position = first_unused_position_after_applying_operation;
        }
        else
        {
            break;
        }
    }

    if (operation_it == operations_.end())
    {
        // Throw error (test first)
    }

    if (first_unused_position == reference_position)
    {
        list<Operation> suffix_operations;
        suffix_operations.splice(suffix_operations.begin(), operations_, operation_it, operations_.end());

        updateCounts();
        return Alignment(first_unused_position, suffix_operations);
    }
    else
    {
        const size_t first_piece_reference_length = reference_position - first_unused_position;
        OperationPair prefix_suffix = splitByReferenceLength(*operation_it, first_piece_reference_length);

        list<Operation> suffix_operations;
        suffix_operations.splice(suffix_operations.begin(), operations_, ++operation_it, operations_.end());
        suffix_operations.push_front(prefix_suffix.second);

        operations_.pop_back();
        operations_.push_back(prefix_suffix.first);

        updateCounts();
        return Alignment(reference_position, suffix_operations);
    }
}

void Alignment::reverse(size_t reference_length)
{
    reference_start_ = reference_length - reference_start_ - referenceLength();
    std::reverse(operations_.begin(), operations_.end());
}

bool Alignment::operator<(const Alignment& other) const
{
    if (reference_start_ != other.reference_start_)
    {
        return reference_start_ < other.reference_start_;
    }

    return operations_ < other.operations_;
}

std::ostream& operator<<(std::ostream& os, const Alignment& alignment)
{
    os << "Ref start: " << alignment.referenceStart() << ", ";
    for (const Operation& operation : alignment)
    {
        os << operation;
    }

    return os;
}
}
