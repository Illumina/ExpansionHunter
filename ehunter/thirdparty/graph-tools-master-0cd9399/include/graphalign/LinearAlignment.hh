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
