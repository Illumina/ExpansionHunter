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

#pragma once

#include <cstdint>
#include <string>

namespace graphtools
{
enum class OperationType
{
    kMatch,
    kMismatch,
    kInsertionToRef,
    kDeletionFromRef,
    kSoftclip,
    kMissingBases
};

// Represents a single alignment operation
class Operation
{
public:
    Operation(OperationType type, uint32_t length)
        : type_(type)
        , length_(length)
    {
    }
    explicit Operation(std::string cigar);

    OperationType type() const { return type_; }
    uint32_t length() const { return length_; }
    uint32_t referenceLength() const;
    uint32_t queryLength() const;

    bool operator==(const Operation& other) const { return type_ == other.type_ && length_ == other.length_; }
    bool operator<(const Operation& other) const;

    std::string generateCigar() const;

private:
    OperationType type_;
    uint32_t length_;
};

OperationType decodeOperationType(char type_encoding);

std::ostream& operator<<(std::ostream& os, OperationType operation_type);
std::ostream& operator<<(std::ostream& os, const Operation& operation);
}
