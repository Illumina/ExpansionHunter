//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "common/ReferenceContigInfo.hh"

namespace ehunter
{

// Represents a contiguous region of a genome using 0-based half-open coordinates
class GenomicRegion
{
public:
    friend std::ostream& operator<<(std::ostream& out, const GenomicRegion& region);

    GenomicRegion(const int32_t contigIndex, int64_t start, int64_t end);

    bool operator<(const GenomicRegion& other) const;

    bool overlaps(const GenomicRegion& other) const;
    int64_t distance(const GenomicRegion& other) const;

    int32_t contigIndex() const { return contigIndex_; }
    int64_t start() const { return start_; }
    int64_t end() const { return end_; }
    int64_t length() const { return end_ - start_; }

    void setContigId(int32_t contigIndex) { contigIndex_ = contigIndex; }
    void setStart(int64_t start) { start_ = start; }
    void setEnd(int64_t end) { end_ = end; }

    bool operator==(const GenomicRegion& other) const
    {
        return contigIndex_ == other.contigIndex_ && start_ == other.start_ && end_ == other.end_;
    }

    bool operator!=(const GenomicRegion& other) const { return !(*this == other); }

    GenomicRegion extend(int length) const;

private:
    int32_t contigIndex_;
    int64_t start_;
    int64_t end_;
};

using GenomicRegionCatalog = std::unordered_map<std::string, GenomicRegion>;

std::ostream& operator<<(std::ostream& out, const GenomicRegion& region);
std::vector<GenomicRegion> merge(std::vector<GenomicRegion> regions, int maxMergeDist = 500);

std::string encode(const ReferenceContigInfo& contigInfo, const GenomicRegion& region);
GenomicRegion decode(const ReferenceContigInfo& contigInfo, const std::string& encoding);

}
