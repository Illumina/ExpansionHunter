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
