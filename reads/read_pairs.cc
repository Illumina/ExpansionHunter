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

#include "reads/read_pairs.h"

#include <stdexcept>

using std::string;
using std::vector;

namespace reads
{

void ReadPairs::Add(const Read& read)
{
    ReadPair& read_pair = read_pairs_[read.fragmentId()];
    const int32_t num_mates_original
        = static_cast<int32_t>(read_pair.first_mate.isSet()) + static_cast<int32_t>(read_pair.second_mate.isSet());

    if (read.is_first_mate && !read_pair.first_mate.isSet())
    {
        read_pair.first_mate = read;
    }

    if (read.isSecondMate() && !read_pair.second_mate.isSet())
    {
        read_pair.second_mate = read;
    }

    const int32_t num_mates_after_add
        = static_cast<int32_t>(read_pair.first_mate.isSet()) + static_cast<int32_t>(read_pair.second_mate.isSet());

    num_reads_ += num_mates_after_add - num_mates_original;
}

const ReadPair& ReadPairs::operator[](const string& fragment_id) const
{
    if (read_pairs_.find(fragment_id) == read_pairs_.end())
    {
        throw std::logic_error("Fragment " + fragment_id + " does not exist");
    }
    return read_pairs_.at(fragment_id);
}

ReadIdToReadReference ReadPairs::GetReads()
{
    ReadIdToReadReference read_references;

    for (auto& kv : read_pairs_)
    {
        ReadPair& read_pair = kv.second;
        if (read_pair.first_mate.isSet())
        {
            const string& read_id = read_pair.first_mate.read_id;
            std::reference_wrapper<Read> read_ref = read_pair.first_mate;
            read_references.emplace(std::make_pair(read_id, read_ref));
        }
        if (read_pair.second_mate.isSet())
        {
            const string& read_id = read_pair.second_mate.read_id;
            std::reference_wrapper<Read> read_ref = read_pair.second_mate;
            read_references.emplace(std::make_pair(read_id, read_ref));
        }
    }

    return read_references;
}

/*
vector<Read::SharedPtr> ReadPairs::GetUnprocessedReads() const
{
    vector<Read::SharedPtr> read_ptrs;

    for (const auto& kv : read_pairs_)
    {
        ReadPair read_pair = kv.second;

        if (read_pair.IsFirstMateSet() && read_pair.FirstMate().CanonicalAlignmentType() == AlignmentType::kUnprocessed)
        {
            read_ptrs.push_back(read_pair.first_mate_ptr_);
        }
        if (read_pair.IsSecondMateSet()
            && read_pair.SecondMate().CanonicalAlignmentType() == AlignmentType::kUnprocessed)
        {
            read_ptrs.push_back(read_pair.second_mate_ptr_);
        }
    }

    return read_ptrs;
} */

void ReadPairs::Clear()
{
    read_pairs_.clear();
    num_reads_ = 0;
}

bool operator==(const ReadPair& read_pair_a, const ReadPair& read_pair_b)
{
    const bool are_first_mates_equal = read_pair_a.first_mate == read_pair_b.first_mate;
    const bool are_second_mates_equal = read_pair_a.second_mate == read_pair_b.second_mate;
    return (are_first_mates_equal && are_second_mates_equal);
}

} // namespace reads
