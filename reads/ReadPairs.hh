//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#pragma once

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include "reads/Read.hh"

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
