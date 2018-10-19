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
