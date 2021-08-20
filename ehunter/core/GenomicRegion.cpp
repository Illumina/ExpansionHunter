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

#include "core/GenomicRegion.hh"

#include <limits>
#include <sstream>
#include <stdexcept>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

using std::istream;
using std::ostream;
using std::string;
using std::unordered_map;
using std::vector;

namespace ehunter
{

GenomicRegion::GenomicRegion(int32_t contigIndex, int64_t start, int64_t end)
    : contigIndex_(contigIndex)
    , start_(start)
    , end_(end)
{
}

bool GenomicRegion::operator<(const GenomicRegion& other) const
{
    if (contigIndex_ != other.contigIndex_)
    {
        return contigIndex_ < other.contigIndex_;
    }

    if (start_ != other.start_)
    {
        return start_ < other.start_;
    }

    return end_ < other.end_;
}

bool GenomicRegion::overlaps(const GenomicRegion& other) const
{
    if (contigIndex_ != other.contigIndex_)
    {
        return false;
    }

    const int64_t leftBound = start_ > other.start_ ? start_ : other.start_;
    const int64_t rightBound = end_ < other.end_ ? end_ : other.end_;

    return leftBound <= rightBound;
}

int64_t GenomicRegion::distance(const GenomicRegion& other) const
{
    if (contigIndex_ != other.contigIndex_)
    {
        return std::numeric_limits<int64_t>::max();
    }

    if (end_ < other.start_)
    {
        return other.start_ - end_;
    }

    if (other.end_ < start_)
    {
        return start_ - other.end_;
    }

    return 0;
}

vector<GenomicRegion> merge(vector<GenomicRegion> regions, int maxMergeDist)
{
    if (regions.empty())
    {
        return regions;
    }

    std::sort(regions.begin(), regions.end());

    GenomicRegion mergedRegion = regions.front();
    vector<GenomicRegion> mergedRegions;

    for (const auto& currentRegion : regions)
    {
        if (currentRegion.distance(mergedRegion) <= maxMergeDist)
        {
            const int64_t furthestEnd = std::max<int64_t>(mergedRegion.end(), currentRegion.end());
            mergedRegion.setEnd(furthestEnd);
        }
        else
        {
            mergedRegions.push_back(mergedRegion);
            mergedRegion = currentRegion;
        }
    }

    if (mergedRegions.empty() || (mergedRegions.back() != mergedRegion))
    {
        mergedRegions.push_back(mergedRegion);
    }

    return mergedRegions;
}

// Returns the range extended by flankSize upstream and downstream.
// NOTE: The right boundary of the extended region may stick past chromosome
// end.
GenomicRegion GenomicRegion::extend(int length) const
{
    const int64_t newStart = start_ >= length ? (start_ - length) : 0;
    const int64_t newEnd = end_ + length;
    return GenomicRegion(contigIndex_, newStart, newEnd);
}

std::ostream& operator<<(std::ostream& out, const GenomicRegion& region)
{
    out << "(" << region.contigIndex_ << "):" << region.start_ << "-" << region.end_;
    return out;
}

string encode(const ReferenceContigInfo& contigInfo, const GenomicRegion& region)
{
    const auto& contigName = contigInfo.getContigName(region.contigIndex());
    return contigName + ":" + std::to_string(region.start()) + "-" + std::to_string(region.end());
}

GenomicRegion decode(const ReferenceContigInfo& contigInfo, const string& encoding)
{
    vector<string> components;
    boost::algorithm::split(components, encoding, boost::algorithm::is_any_of(":-"));

    if (components.size() != 3)
    {
        throw std::logic_error("Unexpected range format: " + encoding);
    }

    const auto& contigName = components[0];
    int32_t contigIndex = contigInfo.getContigId(contigName);

    int64_t start = std::stoi(components[1]);
    int64_t end = std::stoi(components[2]);

    return GenomicRegion(contigIndex, start, end);
}

}
