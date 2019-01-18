//
// Expansion Hunter
// Copyright (c) 2018 Illumina, Inc.
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
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ehunter
{

// Handles translation between contig names and indexes
class ReferenceContigInfo
{
public:
    explicit ReferenceContigInfo(std::vector<std::pair<std::string, int64_t>> namesAndSizes);

    int32_t numContigs() const { return namesAndSizes_.size(); }
    const std::string& getContigName(int32_t contigIndex) const;
    int64_t getContigSize(int32_t contigIndex) const;
    int32_t getContigId(const std::string& contigName) const;

private:
    void assertValidIndex(int32_t contigIndex) const;

    std::vector<std::pair<std::string, int64_t>> namesAndSizes_;
    std::unordered_map<std::string, int32_t> nameToIndex_;
};

std::ostream& operator<<(std::ostream& out, const ReferenceContigInfo& contigInfo);

}