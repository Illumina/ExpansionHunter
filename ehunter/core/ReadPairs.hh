//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
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
//

#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include "core/Read.hh"

namespace ehunter
{
using ReadIdToReadReference = std::unordered_map<std::string, std::reference_wrapper<Read>>;

struct ReadPair
{
    int numMatesSet() const
    {
        return static_cast<int>(firstMate != boost::none) + static_cast<int>(secondMate != boost::none);
    }

    boost::optional<Read> firstMate;
    boost::optional<Read> secondMate;
};

bool operator==(const ReadPair& readPair_a, const ReadPair& readPair_b);

/**
 * Read pair container class
 */
class ReadPairs
{
public:
    typedef std::unordered_map<std::string, ReadPair>::const_iterator const_iterator;
    typedef std::unordered_map<std::string, ReadPair>::iterator iterator;
    const_iterator begin() const { return readPairs_.begin(); }
    const_iterator end() const { return readPairs_.end(); }
    iterator begin() { return readPairs_.begin(); }
    iterator end() { return readPairs_.end(); }

    ReadPairs() = default;
    void Clear();
    void Add(Read read);
    void AddMateToExistingRead(Read mate);

    const ReadPair& operator[](const std::string& fragmentId) const;

    int32_t NumReads() const { return numReads_; }
    int32_t NumCompletePairs() const;

    bool operator==(const ReadPairs& other) const
    {
        return (readPairs_ == other.readPairs_ && numReads_ == other.numReads_);
    }

private:
    std::unordered_map<std::string, ReadPair> readPairs_;
    int32_t numReads_ = 0;
};

}
