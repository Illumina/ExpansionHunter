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

#include "reads/read.h"

namespace reads
{

using ReadIdToReadReference = std::unordered_map<std::string, std::reference_wrapper<Read>>;

struct ReadPair
{
    Read first_mate;
    Read second_mate;
};

bool operator==(const ReadPair& read_pair_a, const ReadPair& read_pair_b);

/**
 * Read pair container class
 */
class ReadPairs
{
public:
    typedef std::unordered_map<std::string, ReadPair>::const_iterator const_iterator;
    typedef std::unordered_map<std::string, ReadPair>::iterator iterator;
    const_iterator begin() const { return read_pairs_.begin(); }
    const_iterator end() const { return read_pairs_.end(); }
    iterator begin() { return read_pairs_.begin(); }
    iterator end() { return read_pairs_.end(); }

    ReadPairs() = default;
    void Clear();
    void Add(const Read& read);

    const ReadPair& operator[](const std::string& fragment_id) const;

    int32_t NumReads() const { return num_reads_; }

    ReadIdToReadReference GetReads();
    // std::vector<Read::SharedPtr> GetUnprocessedReads() const;

    bool operator==(const ReadPairs& other) const
    {
        return (read_pairs_ == other.read_pairs_ && num_reads_ == other.num_reads_);
    }

private:
    std::unordered_map<std::string, ReadPair> read_pairs_;
    int32_t num_reads_ = 0;
};

} // namespace reads