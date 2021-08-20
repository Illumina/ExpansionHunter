//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/ReferenceContigInfo.hh"

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
