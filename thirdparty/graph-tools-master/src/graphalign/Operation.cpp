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

#include "graphalign/Operation.hh"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>

using std::logic_error;
using std::map;
using std::string;
using std::to_string;

namespace graphtools
{

string Operation::generateCigar() const
{
    string cigar_string = to_string(length_);

    std::ostringstream os;
    os << type_;
    cigar_string += os.str();

    return cigar_string;
}

OperationType decodeOperationType(char type_encoding)
{
    switch (type_encoding)
    {
    case 'M':
        return OperationType::kMatch;
    case 'N':
        return OperationType::kMissingBases;
    case 'X':
        return OperationType::kMismatch;
    case 'I':
        return OperationType::kInsertionToRef;
    case 'D':
        return OperationType::kDeletionFromRef;
    case 'S':
        return OperationType::kSoftclip;
    default:
        throw logic_error(to_string(type_encoding) + " is unknown CIGAR operation");
    }
}

Operation::Operation(string cigar)
{
    type_ = decodeOperationType(cigar.back());
    cigar.pop_back();
    length_ = std::stoi(cigar);
}

uint32_t Operation::referenceLength() const
{
    switch (type_)
    {
    case OperationType::kMatch:
    case OperationType::kMismatch:
    case OperationType::kMissingBases:
    case OperationType::kDeletionFromRef:
        return length_;
    default:
        return 0;
    }
}

uint32_t Operation::queryLength() const
{
    switch (type_)
    {
    case OperationType::kMatch:
    case OperationType::kMismatch:
    case OperationType::kMissingBases:
    case OperationType::kInsertionToRef:
    case OperationType::kSoftclip:
        return length_;
    default:
        return 0;
    }
}

bool Operation::operator<(const Operation& other) const
{
    if (type_ != other.type_)
    {
        return type_ < other.type_;
    }

    return length_ < other.length_;
}

std::ostream& operator<<(std::ostream& os, OperationType operation_type)
{
    switch (operation_type)
    {
    case OperationType::kMatch:
        os << 'M';
        break;
    case OperationType::kMismatch:
        os << 'X';
        break;
    case OperationType::kInsertionToRef:
        os << 'I';
        break;
    case OperationType::kDeletionFromRef:
        os << 'D';
        break;
    case OperationType::kSoftclip:
        os << 'S';
        break;
    case OperationType::kMissingBases:
        os << 'N';
    }

    return os;
}

std::ostream& operator<<(std::ostream& os, const Operation& operation)
{
    os << operation.length() << operation.type();
    return os;
}
}
