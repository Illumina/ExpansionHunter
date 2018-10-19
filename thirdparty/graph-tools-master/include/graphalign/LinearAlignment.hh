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

#include <list>
#include <string>

#include "graphalign/Operation.hh"

namespace graphtools
{

// Represents a linear alignment
class Alignment
{
public:
    using size_type = size_t;
    using const_iterator = std::list<Operation>::const_iterator;

    Alignment(int32_t reference_start, std::list<Operation> operations)
        : reference_start_(reference_start)
        , operations_(std::move(operations))
    {
        updateCounts();
    }
    Alignment(uint32_t reference_start, const std::string& cigar);
    const std::list<Operation>& operations() const { return operations_; }
    size_type numOperations() const { return operations_.size(); }
    uint32_t queryLength() const;
    uint32_t referenceLength() const;
    uint32_t referenceStart() const { return reference_start_; }
    void setReferenceStart(uint32_t reference_start) { reference_start_ = reference_start; }

    const_iterator begin() const { return operations_.begin(); }
    const_iterator end() const { return operations_.end(); }
    const Operation& front() const { return operations_.front(); }
    const Operation& back() const { return operations_.back(); }
    size_type size() const { return operations_.size(); }
    bool operator==(const Alignment& other) const
    {
        return operations_ == other.operations_ && reference_start_ == other.reference_start_;
    }
    bool operator<(const Alignment& other) const;

    size_t numMatched() const { return matched_; }
    size_t numMismatched() const { return mismatched_; }
    size_t numClipped() const { return clipped_; }
    size_t numInserted() const { return inserted_; }
    size_t numDeleted() const { return deleted_; }
    std::string generateCigar() const;
    /**
     * Reverses an alignment
     *
     * @param reference_length: Total length of the reference sequence
     */
    void reverse(size_t reference_length);

    /**
     * Splits off a piece of the alignment at the given reference position
     *
     * @param reference_position: Position at which the alignment is to be split
     * @return Suffix alignment
     */
    Alignment splitAtReferencePosition(size_t reference_position);

protected:
    void decodeCigar(const std::string& encoding);
    void updateCounts();

private:
    size_t matched_ = 0;
    size_t mismatched_ = 0;
    size_t clipped_ = 0;
    size_t inserted_ = 0;
    size_t deleted_ = 0;
    size_t missing_ = 0;
    int32_t reference_start_ = 0;
    std::list<Operation> operations_;
};

std::ostream& operator<<(std::ostream& os, const Alignment& alignment);
}
